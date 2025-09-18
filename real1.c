#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <portaudio.h>
#include <sndfile.h>
#include <fftw3.h>
#include <unistd.h> // For sleep()

#define SAMPLE_RATE 44100
#define NUM_CHANNELS 1
#define FRAMES_PER_BUFFER 1024
#define RECORD_SECONDS 5
#define NUM_CHANNELS 16
#define LOW_FREQ 200
#define HIGH_FREQ 7000
#define PULSE_RATE 900

typedef struct {
    float *samples;
    int index;
    int max_samples;
} AudioData;

/* Callback function for audio recording */
static int recordCallback(const void *inputBuffer, void *outputBuffer,
                          unsigned long framesPerBuffer,
                          const PaStreamCallbackTimeInfo *timeInfo,
                          PaStreamCallbackFlags statusFlags, void *userData) {
    AudioData *data = (AudioData *)userData;
    const float *input = (const float *)inputBuffer;

    if (inputBuffer == NULL) return paContinue;

    int samplesToCopy = framesPerBuffer * NUM_CHANNELS;
    if (data->index + samplesToCopy > data->max_samples) return paComplete;

    for (int i = 0; i < samplesToCopy; i++) {
        data->samples[data->index++] = input[i];
    }
    return paContinue;
}

/* Record audio function */
void recordAudio(float *buffer, int sampleRate, int duration) {
    Pa_Initialize();
    PaStream *stream;
    AudioData data = {buffer, 0, sampleRate * duration};

    Pa_OpenDefaultStream(&stream, NUM_CHANNELS, 0, paFloat32, sampleRate,
                         FRAMES_PER_BUFFER, recordCallback, &data);
    Pa_StartStream(stream);
    printf("Recording...\n");
    Pa_Sleep(duration * 1000);
    Pa_StopStream(stream);
    Pa_CloseStream(stream);
    Pa_Terminate();
    printf("Recording complete.\n");
}

/* Play back audio */
void playAudio(float *buffer, int sampleRate, int numSamples) {
    Pa_Initialize();
    PaStream *stream;

    Pa_OpenDefaultStream(&stream, 0, NUM_CHANNELS, paFloat32, sampleRate,
                         FRAMES_PER_BUFFER, NULL, NULL);
    Pa_StartStream(stream);
    Pa_WriteStream(stream, buffer, numSamples);
    Pa_StopStream(stream);
    Pa_CloseStream(stream);
    Pa_Terminate();
}

/* Normalize audio */
void normalizeAudio(float *audio, int length) {
    float maxVal = 0;
    for (int i = 0; i < length; i++) {
        if (fabs(audio[i]) > maxVal) {
            maxVal = fabs(audio[i]);
        }
    }
    if (maxVal > 0) {
        for (int i = 0; i < length; i++) {
            audio[i] /= maxVal;
        }
    }
}

/* Generate biphasic pulse train */
void generatePulseTrain(float *envelope, float *pulseTrain, int length, int sampleRate) {
    int pulseInterval = sampleRate / PULSE_RATE;
    float pulseShape[2] = {1, -1};

    for (int i = 0; i < length; i++) {
        pulseTrain[i] = 0;
    }

    for (int pos = 0; pos < length; pos += pulseInterval) {
        if (pos + 1 < length) {
            pulseTrain[pos] = envelope[pos] * pulseShape[0];
            pulseTrain[pos + 1] = envelope[pos] * pulseShape[1];
        }
    }
}

/* Apply FFT */
void computeFFT(float *signal, int length, int sampleRate) {
    fftwf_complex *out;
    fftwf_plan p;
    float *in = (float *)fftwf_malloc(sizeof(float) * length);
    out = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * (length / 2 + 1));
    p = fftwf_plan_dft_r2c_1d(length, in, out, FFTW_ESTIMATE);

    for (int i = 0; i < length; i++) {
        in[i] = signal[i];
    }

    fftwf_execute(p);
    
    printf("FFT Magnitude Spectrum:\n");
    for (int i = 0; i < length / 2 + 1; i++) {
        printf("Freq: %d Hz, Magnitude: %f\n", (i * sampleRate) / length, sqrt(out[i][0] * out[i][0] + out[i][1] * out[i][1]));
    }

    fftwf_destroy_plan(p);
    fftwf_free(in);
    fftwf_free(out);
}

/* Main function */
int main() {
    int totalSamples = SAMPLE_RATE * RECORD_SECONDS;
    float *audioBuffer = (float *)malloc(sizeof(float) * totalSamples);
    float *pulseTrain = (float *)malloc(sizeof(float) * totalSamples);

    if (!audioBuffer || !pulseTrain) {
        printf("Memory allocation failed.\n");
        return 1;
    }

    // Step 1: Record audio
    recordAudio(audioBuffer, SAMPLE_RATE, RECORD_SECONDS);

    // Step 2: Play recorded audio
    printf("Playing recorded audio...\n");
    playAudio(audioBuffer, SAMPLE_RATE, totalSamples);
    
    // Step 3: Wait 2 seconds
    sleep(2);

    // Step 4: Normalize audio
    normalizeAudio(audioBuffer, totalSamples);

    // Step 5: Generate a pulse train from the envelope
    generatePulseTrain(audioBuffer, pulseTrain, totalSamples, SAMPLE_RATE);

    // Step 6: Play processed audio
    printf("Playing processed audio...\n");
    playAudio(pulseTrain, SAMPLE_RATE, totalSamples);

    // Step 7: Perform FFT analysis
    computeFFT(pulseTrain, totalSamples, SAMPLE_RATE);

    // Free memory
    free(audioBuffer);
    free(pulseTrain);

    return 0;
}
