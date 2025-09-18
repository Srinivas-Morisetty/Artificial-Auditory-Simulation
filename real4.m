clear; clc;

%% Step 1: Record User's Voice
fs = 44100;  % Sampling rate
bits = 16;   % Bit depth
channels = 1; % Mono recording
recObj = audiorecorder(fs, bits, channels);

disp('Start speaking...');
recordblocking(recObj, 5); % Record for 5 seconds
disp('Recording stopped.');

% Retrieve recorded audio
audio = getaudiodata(recObj);

% Play the recorded audio first
disp('Playing recorded audio...');
sound(audio, fs);
pause(length(audio) / fs + 2);

% Silence detection: Trim silent parts
threshold = 0.01;
start_idx = find(abs(audio) > threshold, 1, 'first');
end_idx = find(abs(audio) > threshold, 1, 'last');

if isempty(start_idx) || isempty(end_idx)
    error('No significant speech detected.');
end

audio = audio(start_idx:end_idx);  % Remove leading and trailing silence

% Normalize audio
max_val = max(abs(audio));
if max_val > 0
    audio = audio / max_val;
else
    error('Audio is silent after trimming.');
end

%% Step 2: Define Filter Bank Parameters
num_channels = 16;
low_freq = 200;
high_freq = 7000;

greenwood = @(f) 1.0 * log10(f / 165.4 + 0.88);
inv_greenwood = @(x) 165.4 * (10.^(x / 1.0) - 0.88);
cf = inv_greenwood(linspace(greenwood(low_freq), greenwood(high_freq), num_channels));

filters = cell(num_channels, 1);
envelopes = zeros(length(audio), num_channels);

for i = 1:num_channels
    % Define band edges
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
    
    % Design bandpass filter
    filters{i} = designfilt('bandpassiir', 'FilterOrder', 4, ...
                            'HalfPowerFrequency1', f1, ...
                            'HalfPowerFrequency2', f2, ...
                            'SampleRate', fs);
                        
    % Apply filtering and extract envelope
    band_signal = filtfilt(filters{i}, audio);
    envelope = abs(hilbert(band_signal));
    envelopes(:, i) = envelope;
end

%% Step 3: Generate Electrical Pulse Trains (Biphasic Pulses)
pulse_rate = 900;
pulse_interval = round(fs / pulse_rate);
pulse_shape = [1, -1];
pulse_length = length(pulse_shape);

electrical_impulses = zeros(length(audio), num_channels);

for ch = 1:num_channels
    env = envelopes(:, ch);
    pulse_positions = 1:pulse_interval:length(audio);

    for pos = pulse_positions
        if pos + pulse_length - 1 <= length(audio)
            amplitude = env(pos);
            electrical_impulses(pos:pos+pulse_length-1, ch) = amplitude * pulse_shape;
        end
    end
end

%% Step 4: Create Composite Electrical Stimulation Signal
composite_electrical = sum(electrical_impulses, 2);

%% Step 5: Play Processed Audio
pause(2);
disp('Playing processed audio...');
sound(composite_electrical, fs);

%% Step 6: FFT Analysis of Signals
NFFT = 2^nextpow2(length(composite_electrical));
f_axis = fs/2 * linspace(0, 1, NFFT/2+1);

% FFT for original audio
Y_audio = fft(audio, NFFT);
P_audio = abs(Y_audio/NFFT);
P_audio = P_audio(1:NFFT/2+1);

% FFT for processed audio
Y_proc = fft(composite_electrical, NFFT);
P_proc = abs(Y_proc/NFFT);
P_proc = P_proc(1:NFFT/2+1);

%% Step 7: Plot Everything in One Figure
figure(1);

% Envelope of Channel 1
subplot(3,1,1);
plot((0:length(audio)-1)/fs, envelopes(:,1));
title('Envelope of Channel 1');
xlabel('Time (s)'); ylabel('Amplitude');

% Electrical Impulses for Channel 1
subplot(3,1,2);
plot(electrical_impulses(:,1));
title('Electrical Impulses for Channel 1');
xlabel('Sample'); ylabel('Amplitude');

% Composite Electrical Stimulation Signal
subplot(3,1,3);
plot(composite_electrical);
title('Composite Electrical Stimulation Signal');
xlabel('Sample'); ylabel('Amplitude');
figure(2);
% Original vs. Processed Audio (Time-Domain)
subplot(3,1,1);
hold on;
plot((0:length(audio)-1) / fs, audio, 'b');
plot((0:length(composite_electrical)-1) / fs, composite_electrical, 'g');
hold off;
title('Original (Blue) vs. Processed (green) Audio');
xlabel('Time (s)');
ylabel('Amplitude');

% Frequency Response of Filters
subplot(3,1,2);
hold on;
for i = 1:num_channels
    [h, f] = freqz(filters{i}, 1024, fs);
    plot(f, abs(h));
end
hold off;
title('Frequency Response of Bandpass Filters');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

% Magnitude Spectrum of Original vs. Processed Audio
subplot(3,1,3);
hold on;
plot(f_axis, P_audio, 'b');
plot(f_axis, P_proc, 'r');
hold off;
title('FFT: Original (Blue) vs. Processed (Red) Audio');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
xlim([0 fs/2]);
%jkkfyhshttdh
figure
plot((0:length(audio)-1) / fs, audio, 'b');
title('Original (Blue)');
xlabel('Time (s)');
ylabel('Amplitude');
disp('Processing complete.');
