# Expt-8: Single Target Radar Tracking using Kalman Filter

## Version-1
```matlab
%======================================================================
% Title : Single Target Radar Tracking using Kalman Filter
% Aim   : Estimate range and velocity of a moving target using noisy
%         radar measurements
% Tool  : MATLAB
%======================================================================

clc; clear; close all;

%% Simulation Parameters
T  = 1;              % Sampling time (s)
N  = 100;            % Number of time steps

%% True Target Motion (Ground Truth)
r0 = 1000;           % Initial range (m)
v0 = 50;             % Initial velocity (m/s)

x_true = zeros(2, N);
x_true(:,1) = [r0; v0];

for k = 2:N
    x_true(1,k) = x_true(1,k-1) + T * x_true(2,k-1);
    x_true(2,k) = v0;
end

%% Radar Measurement Noise
sigma_r = 50;        % Range noise std (m)
sigma_v = 5;         % Velocity noise std (m/s)

R = [sigma_r^2  0;
     0           sigma_v^2];

%% Generate Noisy Measurements
z = x_true + [sigma_r * randn(1,N);
              sigma_v * randn(1,N)];

%% Kalman Filter Initialization
x_est = zeros(2, N);
x_est(:,1) = [900; 40];    % Initial estimate

P = diag([500^2, 20^2]);   % Initial covariance

%% State Transition Model
F = [1 T;
     0 1];

%% Measurement Model
H = eye(2);

%% Process Noise
sigma_a = 2;  % Acceleration noise std
Q = sigma_a^2 * [T^4/4  T^3/2;
                 T^3/2  T^2];

%% Kalman Filter Loop
for k = 2:N

    % Prediction
    x_pred = F * x_est(:,k-1);
    P_pred = F * P * F' + Q;

    % Kalman Gain
    K = P_pred * H' / (H * P_pred * H' + R);

    % Update
    x_est(:,k) = x_pred + K * (z(:,k) - H * x_pred);
    P = (eye(2) - K * H) * P_pred;
end

%% Plot Results

time = 1:N;

figure;
subplot(2,1,1);
plot(time, x_true(1,:), 'k-', 'LineWidth', 2); hold on;
plot(time, z(1,:), 'r.', 'MarkerSize', 10);
plot(time, x_est(1,:), 'b--', 'LineWidth', 2);
xlabel('Time Step');
ylabel('Range (m)');
legend('True Range', 'Measured Range', 'Estimated Range');
title('Radar Range Tracking using Kalman Filter');
grid on;

subplot(2,1,2);
plot(time, x_true(2,:), 'k-', 'LineWidth', 2); hold on;
plot(time, z(2,:), 'r.', 'MarkerSize', 10);
plot(time, x_est(2,:), 'b--', 'LineWidth', 2);
xlabel('Time Step');
ylabel('Velocity (m/s)');
legend('True Velocity', 'Measured Velocity', 'Estimated Velocity');
title('Radar Velocity Tracking using Kalman Filter');
grid on;
```
