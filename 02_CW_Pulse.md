# Continuous Wave Radar and Pulse Radar Systems
- `Aim`: To model a CW Radar to estimate the Doppler Shift and calculate target velocity and to implement pulse radar system to estimate the target range.

## Version 1:
```MATLAB
fc = 10e9; % Carrier frequency (Hz)
c = 3e8; % Speed of light (m/s)
v_target = 100; % Target Velocity (m/s)
lambda = c/fc; % Wavelength (m)
fd = 2 * v_target / lambda; % Doppler Shift (Hz)
T = 1e-3; % Observation time (s)
fs = 1e6; % Sampling Frequency (Hz)
t = 0:(1/fs):T; % Time Vector

% Tx and Rx signals
tx = cos(2*pi*fc*t);
rx = cos(2*pi*(fc+fd)*t);

% Mixing (homodyne detection)
mixed = tx .* rx;

% FFT to estimate Doppler Shift
N = length(mixed);
f = fx * (0:N-1) / N;
spectrum = abs(fft(mixed));

% Plot
figure;
plot(f, spectrum);
xlabel("Frequency(Hz)");
ylabel("Amplitude(dB)");
title("Doppler Spectrum");
xlim([0 2*fd]);

[~, idx] = max(spectrum);
estimated_fc = f(idx);
estimated_velocity = (estimated_fc * lambda) / 2;
fprintf('Estimated Velocity: %.2fm/s', estimated_velocity);
```

## Version-2
```matlab
clc; clear; close all;

% ----- Common Parameters -----
c = 3e8; fc = 10e9; lambda = c/fc; % Speed of light, carrier, wavelength

%% ===== CW Radar: Estimate Velocity =====
v_app_true = 60; % Approaching target velocity (m/s)
v_rec_true = -40; % Receding target velocity (m/s)
fs_cw = 1e6; T_cw = 0.01;
t_cw = 0:1/fs_cw:T_cw - 1/fs_cw;

% Doppler shifts for both cases
fd_app = 2 * v_app_true / lambda; % Positive Doppler (approaching)
fd_rec = 2 * v_rec_true / lambda; % Negative Doppler (receding)

% --- Simulation signals ---
tx_cw = sin(2*pi*fc*t_cw);
rx_cw_app = sin(2*pi*(fc + fd_app)*t_cw);
rx_cw_rec = sin(2*pi*(fc + fd_rec)*t_cw);

% --- FFT-based Doppler estimation for approaching ---
mix_app = tx_cw .* rx_cw_app;
N = 2^nextpow2(length(mix_app));
f = fs_cw * (-N/2:N/2-1)/N;
X_app = fftshift(abs(fft(mix_app, N)));
[~, idx_app] = max(X_app);
fd_est_app = abs(f(idx_app));
v_est_app = fd_est_app * lambda / 2;

% --- FFT-based Doppler estimation for receding ---
mix_rec = tx_cw .* rx_cw_rec;
X_rec = fftshift(abs(fft(mix_rec, N)));
[~, idx_rec] = max(X_rec);
fd_est_rec = abs(f(idx_rec));
v_est_rec = fd_est_rec * lambda / 2;

fprintf('--- CW RADAR DOPPLER ESTIMATION ---\n');
fprintf('Approaching Target:\n');
fprintf(' Doppler Shift = %+0.2f Hz, True Velocity = %+0.2f m/s, Estimated Velocity = %+0.2f m/s\n', ...
        fd_app, v_app_true, v_est_app);
fprintf('Receding Target:\n');
fprintf(' Doppler Shift = %+0.2f Hz, True Velocity = %+0.2f m/s, Estimated Velocity = %+0.2f m/s\n\n', ...
        fd_rec, v_rec_true, -v_est_rec);

%% ===== Pulse Radar: Estimate Range =====
r_true = 1500; fs = 10e6; PRF = 10e3; PRI = 1/PRF; Tp = 10e-6; np = 5;
t = 0:1/fs:np*PRI - 1/fs;
tx = zeros(size(t));

pulse_start = round((0:np-1) * PRI * fs) + 1;
pulse_len = round(Tp * fs);
pulse_idx = 0:pulse_len-1;
all_idx = pulse_start' + pulse_idx; all_idx = all_idx(:);
tx(all_idx) = 1;

delay = round(2*r_true/c * fs);
rx = [zeros(1,delay), tx]; rx = rx(1:length(t));

[~, lag] = max(xcorr(rx, tx));
lag = lag - length(tx);
r_est = lag / fs * c / 2;

fprintf('--- PULSE RADAR RANGE ESTIMATION ---\n');
fprintf('Estimated Range: %.2f meters (True: %.2f meters)\n\n', r_est, r_true);

%% ===== CW Visualization with Lower Frequency (for clarity) =====
fc_vis = 1000; % Visualization frequency (1 kHz)
fd_vis_app = 80; % Visible Doppler for approaching
fd_vis_rec = -50; % Visible Doppler for receding

t_vis = 0:1e-5:2e-2;
tx_vis = sin(2*pi*fc_vis*t_vis);
rx_vis_app = sin(2*pi*(fc_vis + fd_vis_app)*t_vis);
rx_vis_rec = sin(2*pi*(fc_vis + fd_vis_rec)*t_vis);

%% ===== Plotting =====
figure('Name','Radar Simulation Results','NumberTitle','off');

% CW transmitted signal
subplot(4,1,1);
plot(t_vis*1e3, tx_vis, 'b');
title('CW RADAR Transmitted Signal (Sine Wave)');
xlabel('Time (ms)'); ylabel('Amplitude');
grid on; xlim([0 10]);

% CW approaching vs receding
subplot(4,1,2);
plot(t_vis*1e3, rx_vis_app, 'r', 'LineWidth', 1.2); hold on;
plot(t_vis*1e3, rx_vis_rec, 'g--', 'LineWidth', 1.2);
legend(sprintf('Approaching (+%d Hz)', fd_vis_app), sprintf('Receding (%d Hz)', fd_vis_rec));
title('CW RADAR: Doppler Shift (Approaching vs Receding)');
xlabel('Time (ms)'); ylabel('Amplitude');
grid on; xlim([0 10]);

% Pulse radar transmitted
subplot(4,1,3);
plot(t*1e6, tx, 'b');
title('Pulse RADAR Transmitted Signal');
xlabel('Time (μs)'); ylabel('Amplitude');
ylim([-0.2 1.2]); grid on;

% Pulse radar received
subplot(4,1,4);
plot(t*1e6, rx, 'r');
title('Pulse RADAR Received Signal (Delayed Echo)');
xlabel('Time (μs)'); ylabel('Amplitude');
ylim([-0.2 1.2]); grid on;

sgtitle('CW and Pulse Radar Simulation Results');
```
