# Expt-9 : Simulation of MIMO Radar for Target Detection

## Version-1
```matlab
%===========================================================
% Program Title: MIMO Radar Simulation for Target Detection
% Aim          : Simulate a MIMO radar system using orthogonal
%                waveforms and estimate target angle
% Software     : MATLAB
%===========================================================

clc; clear; close all;

%% Radar Parameters
c = 3e8;               % Speed of light (m/s)
fc = 10e9;             % Carrier frequency (Hz)
lambda = c/fc;         % Wavelength
d = lambda/2;          % Antenna spacing

Nt = 4;                % Number of transmit antennas
Nr = 4;                % Number of receive antennas

theta_target = 20;     % Target angle (degrees)
SNR = 15;              % Signal-to-noise ratio (dB)

%% Signal Parameters
fs = 10e6;             % Sampling frequency
T = 20e-6;             % Pulse duration
t = 0:1/fs:T-1/fs;
N = length(t);

%% Generate Orthogonal Transmit Signals (Frequency Orthogonal)
txSignals = zeros(Nt, N);

for i = 1:Nt
    txSignals(i,:) = exp(1j*2*pi*(i*1e5)*t);
end

%% Target Steering Vectors
a_tx = exp(1j*2*pi*d*(0:Nt-1).' * sind(theta_target)/lambda);
a_rx = exp(1j*2*pi*d*(0:Nr-1).' * sind(theta_target)/lambda);

%% Received Signal Formation
rxSignals = zeros(Nr, N);

for r = 1:Nr
    for txx = 1:Nt
        rxSignals(r,:) = rxSignals(r,:) + ...
            a_rx(r) * a_tx(txx) * txSignals(txx,:);
    end
end

%% Add Noise
rxSignals = awgn(rxSignals, SNR, 'measured');

%% Matched Filtering
mfOutput = zeros(Nt*Nr, N);

idx = 1;
for r = 1:Nr
    for txx = 1:Nt
        mfOutput(idx,:) = conv(rxSignals(r,:), ...
            conj(fliplr(txSignals(txx,:))), 'same');
        idx = idx + 1;
    end
end

%% Angle Estimation (Beamforming)
theta_scan = -90:0.1:90;
P = zeros(size(theta_scan));

for k = 1:length(theta_scan)
    theta = theta_scan(k);
    a_virtual = exp(1j*2*pi*d*(0:Nt*Nr-1).' * sind(theta)/lambda);
    R = mfOutput * mfOutput';
    P(k) = real(a_virtual' * R * a_virtual);
end

%% Normalize and Plot
P = P / max(P);

figure;
plot(theta_scan, 10*log10(P),'LineWidth',2);
xlabel('Angle (Degrees)');
ylabel('Normalized Power (dB)');
title('MIMO Radar Angle Spectrum');
grid on;
```
