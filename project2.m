

clear;
clc;


[audio, fs] = audioread('sample.wav');
if ~all(isfinite(audio))
    error('Input audio contains non-finite values. Please check your file.');
end
audio = audio(:,1);  % Use first channel if stereo

% Prevent division by zero during normalization
max_val = max(abs(audio));
if max_val == 0
    error('Audio signal is silent (all zeros). Cannot normalize.');
end
audio = audio / max_val;  % Normalize audio

%% Define Filter Bank Parametersand sub bands
num_channels = 16;      
low_freq = 200;        
high_freq = 7000;      


greenwood = @(f) 1.0 * log10(f / 165.4 + 0.88);
inv_greenwood = @(x) 165.4 * (10.^(x / 1.0) - 0.88);
cf = inv_greenwood(linspace(greenwood(low_freq), greenwood(high_freq), num_channels));


filters = cell(num_channels, 1);
envelopes = zeros(length(audio), num_channels);

for i = 1:num_channels
    % Define band edges based on center frequencies
    if i == 1
        f1 = low_freq;
    else
        f1 = (cf(i-1) + cf(i)) / 2;
    end
    if i == num_channels
        f2 = high_freq;
    else
        f2 = (cf(i) + cf(i+1)) / 2;
    end
    
    % Design a bandpass filter for this channel
    filters{i} = designfilt('bandpassiir', 'FilterOrder', 4, ...
                            'HalfPowerFrequency1', f1, ...
                            'HalfPowerFrequency2', f2, ...
                            'SampleRate', fs);
                        
    % Filter the audio and extract the envelope using Hilbert transform
    band_signal = filtfilt(filters{i}, audio);
    envelope = abs(hilbert(band_signal));
    envelopes(:, i) = envelope;
end

%% Step 4: Generate Electrical Pulse Trains (Biphasic Pulses)

pulse_rate = 900;                   % Pulse rate in pulses per second
pulse_interval = round(fs / pulse_rate); % Samples between pulses
% Define a simple biphasic pulse shape (for simulation purposes)
pulse_shape = [1, -1];              % Biphasic: positive phase then negative phase
pulse_length = length(pulse_shape); % Pulse duration in samples

% Preallocate matrix for electrical impulses (one column per channel)
electrical_impulses = zeros(length(audio), num_channels);

% For each channel, generate a pulse train modulated by the envelope
for ch = 1:num_channels
    env = envelopes(:, ch);         % Envelope for current channel
    pulse_positions = 1:pulse_interval:length(audio);  % Pulse start indices
    
    for pos = pulse_positions
        % Check if pulse fits in the signal length
        if pos + pulse_length - 1 <= length(audio)
            % Use the envelope at the pulse time as the amplitude
            amplitude = env(pos);
            electrical_impulses(pos:pos+pulse_length-1, ch) = amplitude * pulse_shape;
        end
    end
end

%%  Creating a Compsite Elec Stimulation Sig
% In actual cochlear implants, pulses are interleaved to avoid channel overlap.
% Here, we sum across channels for a simple composite representation.
composite_electrical = sum(electrical_impulses, 2);

%% Plot and Listen to the Results

% Plot the envelope for one channel and its electrical pulse train
figure;
subplot(3,1,1);
plot((0:length(audio)-1)/fs, envelopes(:,1));
title('Envelope of Channel 1');
xlabel('Time (s)'); ylabel('Amplitude');

subplot(3,1,2);
plot(electrical_impulses(:,1));
title('Electrical Impulses for Channel 1');
xlabel('Sample'); ylabel('Amplitude');

subplot(3,1,3);
plot(composite_electrical);
title('Composite Electrical Stimulation Signal');
xlabel('Sample'); ylabel('Amplitude');

%  Listen to the composite electrical stimulation signal
%The composite electrical signal is a simulation of pulse patterns
% and may not sound like natural speech.
sound(composite_electrical, fs);

%%  FFT Plot for Composite Signal
NFFT = 2^nextpow2(length(composite_electrical));
f_axis = fs/2*linspace(0,1,NFFT/2+1);
Y_comp = fft(composite_electrical, NFFT);
P_comp = abs(Y_comp/NFFT);
P_comp = P_comp(1:NFFT/2+1);

figure;
plot(f_axis, P_comp, 'm');
title('Magnitude Spectrum of Composite Electrical Signal');
xlabel('Frequency (Hz)'); ylabel('Magnitude');
xlim([0 fs/2]);
