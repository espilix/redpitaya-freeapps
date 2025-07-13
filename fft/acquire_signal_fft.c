#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <complex.h>

#include "rp.h"
#include "kiss_fftr.h"

#define BUFFER_SIZE 16384
#define FFT_SIZE (BUFFER_SIZE/2)
#define SAMPLE_RATE 125000000.0  // 125 MHz base rate

typedef struct {
    double frequency;
    double magnitude;
    double phase;
} fft_bin_t;

// Apply Hanning window to reduce spectral leakage
void apply_hanning_window(float *data, int size) {
    for (int i = 0; i < size; i++) {
        double window = 0.5 * (1 - cos(2 * M_PI * i / (size - 1)));
        data[i] *= window;
    }
}

// Perform FFT and calculate magnitude spectrum
int perform_fft_analysis(float *input_data, fft_bin_t *fft_results, int size, double sample_rate) {
    // Allocate memory for FFT
    kiss_fft_scalar *input = (kiss_fft_scalar*)malloc(size * sizeof(kiss_fft_scalar));
    kiss_fft_cpx *output = (kiss_fft_cpx*)malloc((size/2 + 1) * sizeof(kiss_fft_cpx));
    
    if (!input || !output) {
        printf("Memory allocation failed\n");
        return -1;
    }
    
    // Copy input data and convert to kiss_fft_scalar
    for (int i = 0; i < size; i++) {
        input[i] = (kiss_fft_scalar)input_data[i];
    }
    
    // Create FFT configuration
    kiss_fftr_cfg cfg = kiss_fftr_alloc(size, 0, NULL, NULL);
    if (!cfg) {
        printf("FFT configuration failed\n");
        free(input);
        free(output);
        return -1;
    }
    
    // Perform FFT
    kiss_fftr(cfg, input, output);
    
    // Calculate magnitude and phase for each frequency bin
    double freq_resolution = sample_rate / size;
    for (int i = 0; i < size/2; i++) {
        double real = output[i].r;
        double imag = output[i].i;
        
        fft_results[i].frequency = i * freq_resolution;
        fft_results[i].magnitude = sqrt(real * real + imag * imag);
        fft_results[i].phase = atan2(imag, real);
        
        // Convert to dB
        fft_results[i].magnitude = 20 * log10(fft_results[i].magnitude + 1e-12);
    }
    
    // Cleanup
    kiss_fftr_free(cfg);
    free(input);
    free(output);
    
    return 0;
}

// Find peaks in the spectrum
void find_peaks(fft_bin_t *fft_data, int size, double threshold_db) {
    printf("\nDetected Peaks (above %.1f dB):\n", threshold_db);
    printf("Frequency (Hz)\tMagnitude (dB)\tPhase (rad)\n");
    printf("============================================\n");
    
    for (int i = 1; i < size - 1; i++) {
        // Simple peak detection: check if current point is higher than neighbors
        if (fft_data[i].magnitude > fft_data[i-1].magnitude && 
            fft_data[i].magnitude > fft_data[i+1].magnitude &&
            fft_data[i].magnitude > threshold_db) {
            
            printf("%.2f\t\t%.2f\t\t%.4f\n", 
                   fft_data[i].frequency, 
                   fft_data[i].magnitude, 
                   fft_data[i].phase);
        }
    }
}

int main(int argc, char **argv) {
    // Initialize Red Pitaya
    if (rp_Init() != RP_OK) {
        fprintf(stderr, "Red Pitaya API init failed!\n");
        return -1;
    }
    
    printf("Red Pitaya FFT Analysis\n");
    printf("======================\n");
    
    // Setup signal generator (for testing)
    rp_GenReset();
    rp_GenFreq(RP_CH_1, 10000.0);  // 10 kHz test signal
    rp_GenAmp(RP_CH_1, 0.5);
    rp_GenWaveform(RP_CH_1, RP_WAVEFORM_SINE);
    rp_GenOutEnable(RP_CH_1);
    
    // Allocate buffers
    uint32_t buff_size = BUFFER_SIZE;
    float *buffer = (float *)malloc(buff_size * sizeof(float));
    fft_bin_t *fft_results = (fft_bin_t *)malloc(FFT_SIZE * sizeof(fft_bin_t));
    
    if (!buffer || !fft_results) {
        printf("Memory allocation failed\n");
        rp_Release();
        return -1;
    }
    
    // Setup acquisition
    rp_AcqReset();
    rp_AcqSetDecimation(RP_DEC_1);  // No decimation = 125 MSps
    rp_AcqSetTriggerLevel(RP_CH_1, 0.0);
    rp_AcqSetTriggerDelay(ADC_BUFFER_SIZE/2.0);
    
    // Acquire data
    bool fillState = false;
    rp_AcqStart();
    
    sleep(1);  // Wait for buffer to fill
    
    rp_AcqSetTriggerSrc(RP_TRIG_SRC_CHA_PE);
    rp_acq_trig_state_t state = RP_TRIG_STATE_TRIGGERED;
    
    // Wait for trigger
    while(1) {
        rp_AcqGetTriggerState(&state);
        if(state == RP_TRIG_STATE_TRIGGERED) {
            break;
        }
    }
    
    // Wait for buffer to fill
    while(!fillState) {
        rp_AcqGetBufferFillState(&fillState);
    }
    
    rp_AcqStop();
    rp_AcqGetOldestDataV(RP_CH_1, &buff_size, buffer);
    
    printf("Data acquired: %d samples\n", buff_size);
    
    // Get decimation factor for correct sample rate calculation
    uint32_t decimation;
    rp_AcqGetDecimationFactor(&decimation);
    double actual_sample_rate = SAMPLE_RATE / decimation;
    
    printf("Sample rate: %.2f MHz\n", actual_sample_rate / 1e6);
    printf("Frequency resolution: %.2f Hz\n", actual_sample_rate / buff_size);
    
    // Apply window function
    apply_hanning_window(buffer, buff_size);
    
    // Perform FFT analysis
    if (perform_fft_analysis(buffer, fft_results, buff_size, actual_sample_rate) == 0) {
        printf("FFT analysis completed successfully\n");
        
        // Find and display peaks
        find_peaks(fft_results, FFT_SIZE, -60.0);  // -60 dB threshold
        
        // Optionally save spectrum to file
        FILE *fp = fopen("spectrum.txt", "w");
        if (fp) {
            fprintf(fp, "Frequency(Hz)\tMagnitude(dB)\tPhase(rad)\n");
            for (int i = 0; i < FFT_SIZE; i++) {
                fprintf(fp, "%.2f\t%.2f\t%.4f\n", 
                        fft_results[i].frequency, 
                        fft_results[i].magnitude, 
                        fft_results[i].phase);
            }
            fclose(fp);
            printf("Spectrum saved to spectrum.txt\n");
        }
    } else {
        printf("FFT analysis failed\n");
    }
    
    // Cleanup
    free(buffer);
    free(fft_results);
    rp_Release();
    
    return 0;
}
