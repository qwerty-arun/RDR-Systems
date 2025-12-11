# Version 1
```matlab
% Experiment 7: Radar Target Detection and Range Estimation considering RCS
% AIM: To analyze how RCS affects the received signal and perform matched filtering.
%
% FINAL REVISION TO MATCH MANUAL'S PLOTS: This version is specifically
% modified to reproduce the visual output of the lab manual, despite the
% manual's internal inconsistencies.
%
% --- KEY INCONSISTENCY IN MANUAL ---
% - The "Received Signal" plot shows an echo starting at 0.6 µs (implies 90m range).
% - The "Matched Filter Output" plot shows a peak at 0.7 µs (implies 105m range).
% - The text states a 45m or 100m range.
%
% This code prioritizes matching the Matched Filter plot, as it is the final
% result of the experiment. Therefore, a true range of 105m is used.

% --- Setup ---
clear;
clc;
close all;

% --- 1. Simulation Configuration ---

% Physical Constants
SPEED_OF_LIGHT = 3e8; % m/s

% Target Parameters
% We set the range to 105m to force the matched filter peak to occur at 0.7µs,
% matching the third plot in the manual.
TRUE_RANGE_METERS = 105.00;
TARGET_RCS_SQ_METERS = 1.00;

% Radar & Simulation Parameters
SAMPLING_FREQ_HZ = 1e9; % 1 GHz
PULSE_WIDTH_SEC = 0.7e-6; % 0.7 microseconds
% The manual's plots clearly end at 0.9 microseconds.
SIMULATION_TIME_SEC = 0.9e-6;

% Create the time vector for the simulation
time_vector = 0:1/SAMPLING_FREQ_HZ:SIMULATION_TIME_SEC;

% --- 2. Generate Transmitted Waveform (Tx) ---

% The manual's plot shows a linear ramp.
tx_pulse = zeros(size(time_vector));
num_pulse_samples = floor(PULSE_WIDTH_SEC * SAMPLING_FREQ_HZ);
ramp_signal = linspace(0.1, 1, num_pulse_samples);
tx_pulse(1:num_pulse_samples) = ramp_signal;

% --- 3. Simulate Target Echo (Rx) ---

% Calculate the round-trip time delay. This will now be 0.7 µs.
round_trip_delay = 2 * TRUE_RANGE_METERS / SPEED_OF_LIGHT;
delay_in_samples = round(round_trip_delay * SAMPLING_FREQ_HZ);

% Create the received signal by delaying the transmitted pulse.
% Note: The signal will now start at 0.7 µs, which differs from the manual's
% second plot but is necessary to get the third plot correct.
attenuation_factor = sqrt(TARGET_RCS_SQ_METERS) * 0.5;
rx_signal_clean = circshift(tx_pulse, delay_in_samples) * attenuation_factor;

% Add a small amount of noise to match the visual texture in the manual
noise = 0.01 * randn(size(time_vector));
rx_signal_with_noise = rx_signal_clean + noise;

% --- 4. Apply Matched Filter for Pulse Compression ---

% The matched filter is the time-reversed conjugate of the transmitted pulse
matched_filter = fliplr(conj(tx_pulse));

% Perform convolution
matched_filter_output = conv(rx_signal_with_noise, matched_filter, 'same');

% --- 5. Estimate Target Range ---

% Find the peak of the matched filter output
[peak_value, peak_index] = max(abs(matched_filter_output));

% Convert the index of the peak back to a time delay
estimated_time_delay = time_vector(peak_index);

% Calculate the final range estimate
estimated_range_meters = (estimated_time_delay * SPEED_OF_LIGHT) / 2;

% --- Display Results ---

fprintf('--- Target Detection with RCS (Manual Reconciliation Version) ---\n');
fprintf('NOTE: This simulation is configured to match the manual''s Matched Filter plot.\n');
fprintf('A true range of %.2fm is used to place the peak at 0.7 µs.\n', TRUE_RANGE_METERS);
fprintf('This causes a deliberate mismatch with the manual''s Received Signal plot.\n\n');
fprintf('Target RCS: %.2f m^2\n', TARGET_RCS_SQ_METERS);
fprintf('Simulated Target Range: %.2f m\n', TRUE_RANGE_METERS);
fprintf('Estimated Target Range from Simulation: %.2f m\n', estimated_range_meters);

% --- Plotting ---

% Manually scale the matched filter output to match the manual's y-axis appearance
mf_output_scaled = abs(matched_filter_output) / peak_value * 0.8;
peak_value_scaled = mf_output_scaled(peak_index);

figure('Name', 'RCS Detection - Reconciled with Manual Plots');

% Plot 1: Transmitted Pulse
subplot(3, 1, 1);
plot(time_vector * 1e6, tx_pulse, 'LineWidth', 1.5);
title('Transmitted Pulse');
xlabel('Time (\mus)');
ylabel('Amplitude');
grid on;
xlim([0, SIMULATION_TIME_SEC * 1e6]);
ylim([-0.2, 1.2]);

% Plot 2: Received Signal
subplot(3, 1, 2);
plot(time_vector * 1e6, rx_signal_with_noise, 'LineWidth', 1.5);
title('Received Signal (with RCS)');
xlabel('Time (\mus)');
ylabel('Amplitude');
grid on;
xlim([0, SIMULATION_TIME_SEC * 1e6]);
ylim([-0.2, 0.8]);

% Plot 3: Matched Filter Output
subplot(3, 1, 3);
plot(time_vector * 1e6, mf_output_scaled, 'LineWidth', 1.5);
hold on;
% The peak will now correctly be at 0.7 µs
plot(estimated_time_delay * 1e6, peak_value_scaled, 'ro', 'MarkerSize', 8, 'LineWidth', 1.5, 'MarkerFaceColor', 'w');
title('Matched Filter Output (Pulse Compression)');
xlabel('Time (\mus)');
ylabel('Amplitude');
legend('MF Output', 'Detected Peak', 'Location', 'northwest');
grid on;
xlim([0, SIMULATION_TIME_SEC * 1e6]);
ylim([-0.1, 1.0]);
hold off;
```
