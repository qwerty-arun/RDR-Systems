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