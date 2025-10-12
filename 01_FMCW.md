# Frequency Modulated Continuous Wave (FMCW)
- `Aim`: To detect the range of target using a FMCW radar system by simulating transmitter, receiver and beat signals in MATLAB.
## Version 1:
```matlab
% FMCW Radar Range Detection (single target)
clear; close all; clc;

% Constants
c = 3e8;                     % speed of light (m/s)

% Radar / waveform parameters (typical automotive LRR/SRR scale)
fc = 77e9;                   % carrier frequency (Hz)
B  = 150e6;                  % sweep bandwidth (Hz)
T  = 40e-6;                  % sweep time (s)
S  = B/T;                    % chirp slope (Hz/s)

% Sampling and frames
fs        = 2*B;             % ADC sampling rate (Hz) >= Nyquist
Ns        = round(T*fs);     % samples per chirp (fast-time)
Nchirp    = 64;              % number of chirps (slow-time)
PRI       = T;               % repetition interval (no idle for simplicity)

% Target parameters
R_true    = 75;              % target range (m)
v_true    = 0;               % target radial velocity (m/s), set 0 for pure range

% Derived delays and Doppler
tau = 2*R_true/c;            % round-trip delay (s)
fd  = 2*v_true*fc/c;         % Doppler shift (Hz), zero if v_true=0

% Time axes
t_fast = (0:Ns-1).' / fs;    % fast-time within a chirp
t_slow = (0:Nchirp-1) * PRI; % slow-time across chirps

% Preallocate data cube: [fast-time x slow-time]
tx = zeros(Ns, Nchirp);
rx = zeros(Ns, Nchirp);

% Generate Tx and Rx for each chirp (linear up-chirp per PRI)
for k = 1:Nchirp
    t = t_fast + t_slow(k);
    % FMCW up-chirp baseband (complex exponential)
    phi_tx = 2*pi*(fc*t + 0.5*S*t.^2);
    tx(:,k) = exp(1j*phi_tx);

    % Received echo: delayed and Doppler shifted
    t_del = t - tau;
    % Only keep valid samples (t_del >= 0). For simplicity, zero otherwise.
    valid = t_del >= 0;
    phi_rx = 2*pi*(fc*t_del + 0.5*S*t_del.^2) + 2*pi*fd*t; % Doppler on carrier
    rx_ch = zeros(Ns,1);
    rx_ch(valid) = 0.1 * exp(1j*phi_rx(valid));            % attenuation factor
    rx(:,k) = rx_ch;
end

% Dechirp (mix): beat signal
beat = tx .* conj(rx);

% Windowing to reduce leakage
w = hamming(Ns);
beat_win = beat .* w;

% 1D FFT across fast-time (range)
Nfft = 2^nextpow2(Ns);
BF = fft(beat_win, Nfft, 1);

% Frequency axis and range mapping
f_axis = (0:Nfft-1).' * (fs/Nfft);      % beat frequency bins (positive half)
R_axis = (c * f_axis) / (2*S);          % meters

% Take magnitude and average across chirps
range_profile = mean(abs(BF), 2);

% Keep only useful positive ranges (up to max range)
maxR = 200;                              % display limit
idx = R_axis <= maxR;

% Plot
figure;
plot(R_axis(idx), 20*log10(range_profile(idx)/max(range_profile(idx))));
grid on; xlabel('Range (m)'); ylabel('Amplitude (dB)');
title('FMCW Range Profile (1D FFT)');

% Peak detection (simple)
[~,iPeak] = max(range_profile(idx));
R_est = R_axis(idx);
fprintf('Estimated range: %.2f m (true: %.2f m)\n', R_est(iPeak), R_true);

% Assume you already have tx and rx signals
beat = tx .* conj(rx);   % dechirp

% Plot a snippet of the beat waveform
figure;
plot(real(beat(:,1)));
xlabel('Sample'); ylabel('Amplitude');
title('Beat Signal (time domain)');

% FFT to see beat frequency
Nfft = 2^nextpow2(length(beat(:,1)));
BF = fft(beat(:,1), Nfft);
f_axis = (0:Nfft-1)*(fs/Nfft);

figure;
plot(f_axis, abs(BF));
xlim([0 2e6]); % zoom to low frequencies
xlabel('Frequency (Hz)'); ylabel('|Beat Spectrum|');
title('Beat Frequency Spectrum');
```

## Version 2:
```matlab
clear all; close all; clc;
 
%% Radar Specifications
c = 3e8;              % Speed of light (m/s)
fc = 77e9;            % Carrier frequency (Hz)
R_max = 200;          % Maximum Range (m)
range_res = 1;        % Range Resolution (m)
target_range = 50;    % Actual Target Distance (m)
 
%% FMCW Waveform Parameters
B = c / (2 * range_res);          % Bandwidth (Hz)
T_chirp = 5.5 * 2 * R_max / c;    % Chirp time (s)
slope = B / T_chirp;              % Chirp slope (Hz/s)
 
%% Sampling Parameters
Nd = 128;                        % Number of chirps (slow time dimension)
Nr = 1024;                       % Number of samples per chirp (fast time)
t = linspace(0, Nd*T_chirp, Nr*Nd);  % Total time for all chirps
 
%% Transmit and Receive Signals
Tx = cos(2*pi*(fc*t + (slope * t.^2)/2));   % Transmitted chirp
td = 2 * target_range / c;                   % Time delay for target echo
Rx = cos(2*pi*(fc*(t - td) + (slope * (t - td).^2)/2));  % Received echo
 
%% Beat Signal Generation (De-chirping)
Mix = Tx .* Rx;   % Element-wise multiplication → beat signal
 
%% Reshape the Mixed Signal
Mix_reshape = reshape(Mix, [Nr, Nd]);   % Arrange into [samples per chirp × number of chirps]
 
%% Range FFT (1D FFT across fast time)
signal_fft = fft(Mix_reshape, Nr);      % FFT along range dimension
signal_fft = abs(signal_fft / Nr);      % Normalize magnitude
signal_fft = signal_fft(1:Nr/2, :);     % Take only positive frequencies
 
%% Plot Range FFT for the First Chirp
figure;
plot(signal_fft(:,1));
title('Range from First Chirp');
xlabel('Range Bins');
ylabel('Amplitude');
grid on;
 
%% Calculate Range Axis
range_axis = linspace(0, R_max, Nr/2);
 
%% Plot Detected Range
figure;
plot(range_axis, signal_fft(:,1));
title('Detected Range');
xlabel('Range (m)');
ylabel('Amplitude');
grid on;
 
%% Optional: Display detected peak
[~, peak_index] = max(signal_fft(:,1));
detected_range = range_axis(peak_index);
disp(['Detected Target Range ≈ ', num2str(detected_range), ' m']);
```