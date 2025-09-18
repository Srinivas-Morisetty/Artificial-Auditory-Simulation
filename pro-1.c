#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sndfile.h>
#include <fftw3.h>

#define NUM_CHANNELS 16
#define LOW_FREQ 200.0
#define HIGH_FREQ 7000.0
#define FILTER_ORDER 4
#define PI 3.141592653589793

// Function to compute Greenwood function
double greenwood(double f) {
    return log10(f / 165.4 + 0.88);
}

double inv_greenwood(double x) {
    return 165.4 * (pow(10.0, x) - 0.88);
}

// Function to generate center frequencies using Greenwood mapping
void compute_center_frequencies(double cf[], int num_channels) {
    double low = greenwood(LOW_FREQ);
    double high = greenwood(HIGH_FREQ);

    for (int i = 0; i < num_channels; i++) {
        double x = low + (high - low) * i / (num_channels - 1);
        cf[i] = inv_greenwood(x);
    }
}

// Simple Bandpass filter (IIR approximation)
void bandpass_filter(double* signal, int length, double fs, double f1, double f2) {
    double alpha1 = exp(-2.0 * PI * f1 / fs);
    double alpha2 = exp(-2.0 * PI * f2 / fs);
    
    for (int i = 1; i < length; i++) {
        signal[i] = alpha1 * signal[i] + (1 - alpha1) * signal[i - 1];
        signal[i] = alpha2 * signal[i] + (1 - alpha2) * signal[i - 1];
    }
}

// Hilbert Transform Approximation (Envelope Extraction)
void hilbert_transform(double* signal, double* envelope, int length) {
    fftw_complex *in, *out;
    fftw_plan plan_forward, plan_inverse;

    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * length);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * length);

    for (int i = 0; i < length; i++) {
        in[i][0] = signal[i];
        in[i][1] = 0.0;
    }

    plan_forward = fftw_plan_dft_1d(length, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan_forward);

    for (int i = 1; i < length / 2; i++) {
        out[i][0] *= 2;
        out[i][1] *= 2;
    }
    for (int i = length / 2 + 1; i < length; i++) {
        out[i][0] = 0;
        out[i][1] = 0;
    }

    plan_inverse = fftw_plan_dft_1d(length, out, in, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan_inverse);

    for (int i = 0; i < length; i++) {
        envelope[i] = sqrt(in[i][0] * in[i][0] + in[i][1] * in[i][1]) / length;
    }

    fftw_destroy_plan(plan_forward);
    fftw_destroy_plan(plan_inverse);
    fftw_free(in);
    fftw_free(out);
}

// Generate cosine carrier wave
void generate_carrier(double* carrier, int length, double fs, double freq) {
    for (int i = 0; i < length; i++) {
        carrier[i] = cos(2 * PI * freq * i / fs);
    }
}

// Normalize audio to prevent clipping
void normalize_audio(double* audio, int length) {
    double max_val = 0.0;
    for (int i = 0; i < length; i++) {
        if (fabs(audio[i]) > max_val) max_val = fabs(audio[i]);
    }
    if (max_val > 0) {
        for (int i = 0; i < length; i++) {
            audio[i] /= max_val;
        }
    }
}

// Main Processing
int main() {
    SNDFILE *infile, *outfile;
    SF_INFO sfinfo;
    
    // Open input audio file
    infile = sf_open("sample_audio.wav", SFM_READ, &sfinfo);
    if (!infile) {
        printf("Error opening input file.\n");
        return 1;
    }

    int fs = sfinfo.samplerate;
    int num_samples = sfinfo.frames;
    double *audio = (double*)malloc(num_samples * sizeof(double));

    // Read the audio file
    sf_read_double(infile, audio, num_samples);
    sf_close(infile);

    // Normalize input audio
    normalize_audio(audio, num_samples);

    // Define filter bank parameters
    double cf[NUM_CHANNELS];
    compute_center_frequencies(cf, NUM_CHANNELS);

    // Process each channel
    double *modulated_signals = (double*)calloc(num_samples, sizeof(double));

    for (int i = 0; i < NUM_CHANNELS; i++) {
        // Define band edges
        double f1 = (i == 0) ? LOW_FREQ : (cf[i - 1] + cf[i]) / 2.0;
        double f2 = (i == NUM_CHANNELS - 1) ? HIGH_FREQ : (cf[i] + cf[i + 1]) / 2.0;

        // Filter the audio signal
        double *band_signal = (double*)malloc(num_samples * sizeof(double));
        for (int j = 0; j < num_samples; j++) {
            band_signal[j] = audio[j];
        }
        bandpass_filter(band_signal, num_samples, fs, f1, f2);

        // Extract envelope
        double *envelope = (double*)malloc(num_samples * sizeof(double));
        hilbert_transform(band_signal, envelope, num_samples);

        // Generate carrier wave
        double *carrier = (double*)malloc(num_samples * sizeof(double));
        generate_carrier(carrier, num_samples, fs, cf[i]);

        // Modulate signal
        for (int j = 0; j < num_samples; j++) {
            modulated_signals[j] += envelope[j] * carrier[j];
        }

        // Free memory for this channel
        free(band_signal);
        free(envelope);
        free(carrier);
    }

    // Normalize processed audio
    normalize_audio(modulated_signals, num_samples);

    // Save processed audio
    sfinfo.channels = 1;
    outfile = sf_open("sample-file-2.wav", SFM_WRITE, &sfinfo);
    sf_write_double(outfile, modulated_signals, num_samples);
    sf_close(outfile);


    printf("Processed audio saved as 'sample-file-2.wav'\n");

    // Free memory
    free(audio);
    free(modulated_signals);

    return 0;
}
