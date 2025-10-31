# To simulate and analyze the spectogram of FMCW radar signals for multiple automotive radar systems in MATLAB
## AIM: The aim is to simulate a pulse radar system and generate the transmitted and received signals including target echoes, clutter, and noise and to apply matched filtering to detect the target and accurately estimate its range

### Version-1
```matlab
%% FMCW Radar Simulation with Range Estimation and Composite 3D Ambiguity Function (dB Scale)
clc; clear; close all;

%% -------------------- RADAR PARAMETERS --------------------
c = 3e8;                % Speed of light (m/s)
fc = 77e9;              % Carrier frequency (Hz)
B = 150e6;              % Sweep bandwidth (Hz)
Tchirp = 5.5e-6;        % Chirp duration (s)
Ns = 1024;              % Samples per chirp
Fs = Ns / Tchirp;       % Sampling frequency (Hz)
lambda = c / fc;        % Wavelength

%% -------------------- TARGET PARAMETERS --------------------
R = [50, 100, 150, 500];     % Target ranges (m)
v = [10, -5, 0, 20];         % Target velocities (m/s)
rcs = [1, 0.8, 0.6, 0.5];    % Radar cross-sections (relative)

%% -------------------- DERIVED PARAMETERS --------------------
slope = B / Tchirp;             % Chirp slope (Hz/s)
t = linspace(0, Tchirp, Ns);    % Fast-time vector

%% -------------------- TRANSMITTED SIGNAL --------------------
Tx = exp(1j * 2 * pi * (fc * t + 0.5 * slope * t.^2));  % FMCW up-chirp

%% -------------------- RECEIVED SIGNAL (MULTIPLE TARGETS) --------------------
Rx = zeros(size(Tx));
for k = 1:length(R)
    tau = 2 * R(k) / c;                % Time delay (s)
    fd = 2 * v(k) / lambda;            % Doppler shift (Hz)
    echo = rcs(k) * exp(1j * 2 * pi * (fc*(t - tau) + 0.5 * slope * (t - tau).^2 + fd * t));
    Rx = Rx + echo;
end

%% -------------------- ADD CLUTTER AND NOISE --------------------
clutter = 0.2 * (randn(size(Tx)) + 1j * randn(size(Tx))); % Clutter noise
SNR_dB = 15;
Rx = awgn(Rx + clutter, SNR_dB, 'measured');

%% -------------------- MIXING (DECHIRPING) --------------------
Mix = Tx .* conj(Rx);

%% -------------------- MATCHED FILTERING --------------------
mf = conj(fliplr(Tx));              % Matched filter (time reversed)
mf_output = conv(Mix, mf, 'same');

%% -------------------- RANGE ESTIMATION --------------------
Nfft = 2^nextpow2(Ns);
range_fft = abs(fft(mf_output, Nfft));
range_axis = (c * Fs / (2 * slope)) * (0:Nfft/2 - 1) / Nfft;

figure;
plot(range_axis, range_fft(1:Nfft/2)/max(range_fft), 'b', 'LineWidth', 1.2);
xlabel('Range (m)');
ylabel('Normalized Amplitude');
title('Range Response after Matched Filtering');
grid on;
hold on;

%% ---- Target Range Detection ----
threshold = 0.4 * max(range_fft);  % detection threshold
[~, locs] = findpeaks(range_fft(1:Nfft/2), 'MinPeakHeight', threshold);
estimated_ranges = range_axis(locs);

% Print detected target ranges
disp('-------------------------------------------------');
disp('Detected Target Ranges (in meters):');
disp(estimated_ranges.');
disp('-------------------------------------------------');

% Mark detected targets on the plot
plot(estimated_ranges, range_fft(locs)/max(range_fft), 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);
text(estimated_ranges + 2, range_fft(locs)/max(range_fft), ...
     strcat(string(round(estimated_ranges,1)), ' m'), 'Color', 'r', 'FontSize', 9);
hold off;

%% -------------------- CUSTOM 3D AMBIGUITY FUNCTION PARAMETERS --------------------
disp('Computing 3D ambiguity functions (this may take a few seconds)...');

fd_max = 5e3;                      % Doppler frequency range (Hz)
num_fd = 150;                     % Doppler samples
num_tau = 150;                    % Time delay samples

tau_vec = linspace(-Tchirp/2, Tchirp/2, num_tau);  % Time delay vector (s)
fd_vec  = linspace(-fd_max, fd_max, num_fd);       % Doppler frequency vector (Hz)

%% -------------------- BASE AMBIGUITY FUNCTION (zero delay, zero Doppler) --------------------
AF_base = zeros(num_fd, num_tau);
for i = 1:num_fd
    fd_shift = fd_vec(i);
    for j = 1:num_tau
        tau = tau_vec(j);
        shifted = exp(1j*2*pi*fd_shift*t) .* interp1(t, Tx, t - tau, 'linear', 0);
        AF_base(i,j) = abs(sum(Tx .* conj(shifted)));
    end
end
AF_base = AF_base / max(AF_base(:));  % Normalize base AF

%% -------------------- COMPOSITE AMBIGUITY FUNCTION --------------------
AF_composite = zeros(num_fd, num_tau);

for k = 1:length(R)
    tau_k = 2 * R(k) / c;         % time delay (s)
    fd_k = 2 * v(k) / lambda;     % Doppler freq (Hz)

    % Find closest indices in tau_vec and fd_vec
    [~, idx_tau] = min(abs(tau_vec - tau_k));
    [~, idx_fd] = min(abs(fd_vec - fd_k));

    % Shift AF_base by target's tau and Doppler indices (circular shift)
    AF_shifted = circshift(AF_base, [idx_fd - floor(num_fd/2), idx_tau - floor(num_tau/2)]);

    % Weight by RCS
    AF_composite = AF_composite + rcs(k) * AF_shifted;
end

% Normalize composite AF and convert to dB scale
AF_composite = AF_composite / max(AF_composite(:));
AF_composite_dB = 20 * log10(AF_composite);
AF_composite_dB(AF_composite_dB < -50) = -50;  % Dynamic range limit

%% -------------------- PLOT COMPOSITE 3D AMBIGUITY FUNCTION --------------------
figure;
mesh(tau_vec, fd_vec/1e3, AF_composite_dB);
xlabel('Time Delay (s)');
ylabel('Doppler Frequency (kHz)');
zlabel('Ambiguity Magnitude (dB)');
title('Composite 3D Ambiguity Function of Multiple Targets');
colormap('jet');
shading interp;
view(45, 30);
colorbar;

%% -------------------- SPECTROGRAM --------------------
figure;
spectrogram(Mix, 128, 120, 512, Fs, 'yaxis');
title('Spectrogram of Received FMCW Signal');
ylabel('Frequency (Hz)');
xlabel('Time (s)');
colormap('turbo');
```
