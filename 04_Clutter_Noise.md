# Expt-4: To estimate the range of a target in the presence of clutter and noise using pulse radar and matched filtering in MATLAB

# AIM: To simulate a pulse radar system and generate the tx and rx signals including target echos, clutter, and noise and to apply matched filtering to detect the target and accurately estimate its range.

## Version-1:
```matlab
% Pulse radar simulation with LFM pulse, clutter, noise and matched filtering
% - Generates tx chirp, simulates target echoes and clutter, adds noise
% - Applies matched filter (time-reversed complex-conjugate of tx)
% - Uses CA-CFAR to detect peaks and estimate range
% - Plots results and prints estimated range vs true range
%
% Author: ChatGPT
% Date: 2025-10-24

clear; close all; clc;

%% Physical + waveform parameters
c = 3e8;                  % speed of light (m/s)
fc = 0;                   % carrier for baseband simulation (0 => complex baseband)
B = 2e6;                  % chirp bandwidth (Hz)
Tp = 20e-6;               % pulse duration (s)
fs = 20e6;                % sampling rate (Hz) (>= 2*B recommended). adjust for resolution
PRI = 1e-3;               % pulse repetition interval (s) - not used heavily here (single pulse)
Npad = 4096;              % total fast-time samples in range window (zero-padding size for conv)

%% Derived values
k = B / Tp;               % chirp rate (Hz/s)
Nt = round(Tp * fs);      % number of samples in pulse
t_pulse = (0:Nt-1) / fs;  % pulse time axis (seconds)
range_res = c / (2*B);    % theoretical range resolution
fprintf('Pulse samples: %d, Range resolution: %.2f m\n', Nt, range_res);

%% Create transmit LFM (complex baseband)
f0 = -B/2; % center LFM at 0 (baseband) -> sweep from -B/2 to B/2
tx = exp(1j*2*pi*(f0 .* t_pulse + 0.5 * k .* t_pulse.^2));  % complex chirp
% apply window (optional) to reduce sidelobes
wnd = hamming(Nt).';
tx = tx .* wnd;

% normalize tx energy to 1 (so echo amplitude scaling behaves nicely)
E_tx = sum(abs(tx).^2);
tx = tx / sqrt(E_tx);
E_tx = sum(abs(tx).^2); % should be 1

%% Range axis and receive buffer
N = Npad;
t_fast = (0:N-1) / fs;                   % fast time axis for buffer
range_axis = (c * t_fast) / 2;           % corresponding ranges (m)

%% Targets and clutter configuration
% Single main target
target_range = 2500;     % meters (true target range)
tau_target = 2 * target_range / c;       % two-way delay (s)
sample_delay_target = round(tau_target * fs);

% amplitude loss (free-space approx ~ 1/r^2). We'll set base amplitude:
amp0 = 1;  % reference amplitude at some reference; we'll scale by 1/R^2
amp_target = amp0 / (target_range^2);    % amplitude of echo

% Add some phase for realism
phase_target = exp(1j*2*pi*rand);

% clutter: many weak scatterers across a range interval
num_clutter = 80;
clutter_ranges = sort( (50 + (range_axis(end)-50).*rand(1,num_clutter)) ); % random ranges
clutter_amplitudes = 0.2 .* (0.5 + rand(1,num_clutter)); % weaker than target
clutter_phases = exp(1j*2*pi*rand(1,num_clutter));

%% Build received (noise-free) baseband signal (single pulse)
rx = zeros(1, N);  % complex baseband samples

% Insert target echo: delayed version of tx
if sample_delay_target + Nt <= N
    rx(sample_delay_target + (1:Nt)) = rx(sample_delay_target + (1:Nt)) ...
        + amp_target * phase_target .* tx;
else
    warning('Target echo falls outside receive buffer. Increase Npad or move target closer.');
end

% Insert clutter echoes (each is a scaled, delayed version of tx)
for ii = 1:num_clutter
    r = clutter_ranges(ii);
    tau_c = 2 * r / c;
    d_samp = round(tau_c * fs);
    amp_c = clutter_amplitudes(ii) / (r^2 + 1);  % smaller dependence on range
    ph = clutter_phases(ii);
    if d_samp + Nt <= N
        rx(d_samp + (1:Nt)) = rx(d_samp + (1:Nt)) + amp_c * ph .* tx;
    end
end

%% Additive white Gaussian noise (complex) for desired SNR
SNR_dB = 0;  % SNR in dB as energy ratio (can change to test)
SNR_linear = 10^(SNR_dB/10);

% We normalized tx energy to 1, and target echo amplitude is amp_target, so the signal energy
% at receiver (for the target) is approx E_sig = amp_target^2 * E_tx = amp_target^2.
E_sig = amp_target^2 * E_tx;  % since E_tx ~ 1

% Set noise variance such that SNR = E_sig / sigma2
sigma2 = E_sig / SNR_linear;
if sigma2 == 0
    noise = zeros(1,N);
else
    noise = sqrt(sigma2/2) * (randn(1,N) + 1j*randn(1,N));
end
rx_noisy = rx + noise;

%% Matched filter
mf = conj(flipud(tx(:))).';   % time-reversed complex conjugate
% zero-pad for clean convolution
y = conv(rx_noisy, mf, 'same');   % matched filter output (complex)
y_mag = abs(y);
y_power = y_mag.^2;

% Range axis for matched filter output is same as t_fast (we used 'same')

%% Simple detection: CA-CFAR
% CA-CFAR parameters
num_train = 20;     % training cells each side
num_guard = 4;      % guard cells each side
Pfa = 1e-4;         % desired false alarm probability

mag2 = y_power;     % detection statistic

% Pre-allocate threshold and detection arrays
threshold = zeros(1,N);
detected = false(1,N);

half_win = num_train + num_guard;
for kidx = (half_win+1):(N-half_win)
    % indices for training cells
    lead_train_idx = (kidx - half_win):(kidx - num_guard - 1);
    tail_train_idx = (kidx + num_guard + 1):(kidx + half_win);
    training_cells = mag2([lead_train_idx, tail_train_idx]);
    noise_level = mean(training_cells); % CA average
    % threshold multiplier from Pfa
    T = length(training_cells); % total training cells
    if T > 0
        alpha = T*(Pfa^(-1/T) - 1);    % CA-CFAR multiplier (approx)
    else
        alpha = 1;
    end
    threshold(kidx) = alpha * noise_level;
    if mag2(kidx) > threshold(kidx)
        detected(kidx) = true;
    end
end

% find peaks among detections (local max)
% find indices which are detected and also local maxima
det_idx = find(detected);
local_peaks = [];
for idx = det_idx
    left = max(idx-1,1);
    right = min(idx+1,N);
    if y_mag(idx) == max(y_mag(left:right))
        local_peaks(end+1) = idx; %#ok<SAGROW>
    end
end

% choose the strongest peak (if any)
if isempty(local_peaks)
    fprintf('No detection by CFAR.\n');
    [~, peak_idx] = max(y_mag);
    fprintf('Reporting strongest peak anyway.\n');
else
    [~, imax] = max(y_mag(local_peaks));
    peak_idx = local_peaks(imax);
end

estimated_range = range_axis(peak_idx);
true_range = target_range;

%% Display results
fprintf('True target range: %.2f m\n', true_range);
fprintf('Estimated target range (peak): %.2f m\n', estimated_range);
range_error = estimated_range - true_range;
fprintf('Range error: %.2f m\n', range_error);

%% Plotting
figure('Position',[100 100 1000 800]);

subplot(4,1,1);
plot(t_pulse*1e6, real(tx));
xlabel('Time (us)'); ylabel('Amplitude (real)');
title('Transmit LFM Pulse (real part)');
grid on;

subplot(4,1,2);
plot(range_axis, real(rx_noisy));
xlim([0 max(range_axis)]);
xlabel('Range (m)'); ylabel('Amplitude (real)');
title('Received signal (noisy) — fast-time (real part)');
grid on;

subplot(4,1,3);
plot(range_axis, 20*log10(y_mag + eps));
xlabel('Range (m)'); ylabel('Matched filter magnitude (dB)');
title('Matched-filter output (range response)');
xlim([0 max(range_axis)]);
grid on;
hold on;
% mark true and estimated ranges
plot(true_range, 20*log10(y_mag(round(tau_target*fs)+1)+eps), 'gv', 'MarkerFaceColor','g', 'DisplayName','True target (nominal)');
plot(estimated_range, 20*log10(y_mag(peak_idx)+eps), 'r^', 'MarkerFaceColor','r', 'DisplayName','Estimated peak');
legend('Matched filter output','True range (nom)','Estimated peak');

subplot(4,1,4);
plot(range_axis, 10*log10(threshold+eps));
hold on;
plot(range_axis, 20*log10(y_mag+eps));
xlabel('Range (m)'); ylabel('dB');
title('Matched-filter output and CFAR threshold');
xlim([0 max(range_axis)]);
legend('CFAR threshold','Matched filter');
grid on;

% Mark detections on the range plot
figure('Position',[200 200 800 300]);
plot(range_axis, 20*log10(y_mag+eps),'LineWidth',1.2);
hold on;
plot(range_axis(detected), 20*log10(y_mag(detected)+eps), 'ro', 'MarkerFaceColor','r');
plot(estimated_range, 20*log10(y_mag(peak_idx)+eps), 'kp', 'MarkerFaceColor','k', 'MarkerSize',10);
xlabel('Range (m)'); ylabel('Magnitude (dB)');
title('Detections (CFAR)');
xlim([0 max(range_axis)]);
legend('Matched filter output','CFAR detections','Selected peak');
grid on;

%% Print a small table of the top detected peaks (range and SNR-like amplitude)
fprintf('\nTop detected peaks (range, amplitude):\n');
% list top 5 local maxima by amplitude
[~, sort_idx] = sort(y_mag, 'descend');
topN = 5;
cnt = 0;
printed = 0;
for s = 1:length(sort_idx)
    idx = sort_idx(s);
    % only print if it's a local maxima
    left = max(idx-1,1);
    right = min(idx+1,N);
    if y_mag(idx) == max(y_mag(left:right))
        printed = printed + 1;
        fprintf('  Peak %d: range = %.2f m, mag = %.4e\n', printed, range_axis(idx), y_mag(idx));
        if printed >= topN, break; end
    end
end

%% End of script
```

## Version-2
```matlab
% Pulse Radar Simulation - Range Estimation via Matched Filtering
clear; close all; clc;

c = 3e8;           % Speed of light (m/s)
fs = 20e6;         % Sampling rate (Hz)
Tp = 20e-6;        % Pulse duration (s)
B = 2e6;           % Bandwidth (Hz)
k = B/Tp;          % Chirp rate
t = 0:1/fs:Tp-1/fs;
tx = exp(1j*pi*k*t.^2);  % LFM pulse (baseband)
tx = tx / sqrt(sum(abs(tx).^2)); % Normalize

% Target and clutter setup
target_range = 2500;                     % meters
tau_t = 2*target_range/c;                % delay
samp_delay = round(tau_t*fs);            % sample delay
N = 4096;
rx = zeros(1,N);
rx(samp_delay+(1:length(tx))) = 1e-5 * exp(1j*2*pi*rand) * tx; % target echo

% Add simple clutter (random echoes)
for i = 1:10
    d = randi([50 2000]);  % random range
    rx(round(2*d/c*fs)+(1:length(tx))) = rx(round(2*d/c*fs)+(1:length(tx))) + 2e-6*exp(1j*2*pi*rand)*tx;
end

% Add noise
rx = rx + (randn(1,N)+1j*randn(1,N))*1e-6;

% Matched filter
y = abs(conv(rx, conj(flip(tx)), 'same'));
range_axis = (0:N-1)*c/(2*fs);

% Find peak
[~, idx] = max(y);
est_range = range_axis(idx);
fprintf('True range: %.1f m | Estimated range: %.1f m\n', target_range, est_range);

% Plot
subplot(3,1,1); plot(real(tx)); title('Transmit Pulse (Real Part)');
subplot(3,1,2); plot(real(rx)); title('Received Signal (Real Part)');
subplot(3,1,3); plot(range_axis, y); hold on;
plot(est_range, y(idx), 'ro'); title('Matched Filter Output'); xlabel('Range (m)');
```

## Version-3
```matlab
clear; clc; close all;

%% ====================== RADAR PARAMETERS ======================
c = 3e8;                   % Speed of light (m/s)
Tp = 1e-6;                 % Pulse width (1 µs)
fs = 10e6;                  % Sampling frequency (10 MHz)
R_target = 15e3;           % Target range (15 km)
SNR_dB = 10;                % Signal-to-noise ratio (dB)

% Multiple clutter sources (in meters)
R_clutter = [5e3, 8e3, 12e3];       % Example: 5 km, 8 km, 12 km
clutter_amplitude = [0.8, 0.6, 0.4]; % Reflection strengths

%% ====================== DERIVED PARAMETERS ======================
N = round(Tp * fs);         % Samples per pulse

%% ====================== TRANSMITTED RECTANGULAR PULSE ======================
% Create a long time vector to clearly see the full pulse
t_tx = (0:5*N-1)/fs;         % 5 times the pulse duration
tx_rect = zeros(1, length(t_tx));
tx_rect(1:N) = 1;             % Pulse starts at t = 0

%% ====================== RECEIVED SIGNAL (TARGET + MULTIPLE CLUTTER + NOISE) ======================
% Convert ranges to delays (seconds)
delay_target  = 2 * R_target / c;
delay_clutter = 2 * R_clutter / c;

% Convert delays to sample numbers
n_delay_target  = round(delay_target  * fs);
n_delay_clutter = round(delay_clutter * fs);

% Create received signal buffer (long enough to contain all echoes)
rx_length = length(t_tx) + n_delay_target + 200;
rx_signal = zeros(1, rx_length);

% Add clutter echoes
for i = 1:length(R_clutter)
    rx_signal(n_delay_clutter(i) + (1:N)) = ...
        rx_signal(n_delay_clutter(i) + (1:N)) + clutter_amplitude(i) * tx_rect(1:N);
end

% Add target echo
rx_signal(n_delay_target + (1:N)) = rx_signal(n_delay_target + (1:N)) + 1 * tx_rect(1:N);

% Add AWGN noise
rx_signal = awgn(rx_signal, SNR_dB, 'measured');

% Time axis for received signal
t_rx = (0:rx_length-1)/fs;   % seconds

%% ====================== MATCHED FILTERING ======================
mf = flip(tx_rect(1:N));             % Matched filter (time-reversed)
mf_output = conv(rx_signal, mf, 'same');
t_mf = (0:length(mf_output)-1)/fs;  % seconds

%% ====================== TARGET & CLUTTER DETECTION ======================
[pks, locs] = findpeaks(abs(mf_output), 'MinPeakHeight', 0.3*max(abs(mf_output)), ...
    'MinPeakDistance', round(0.5e-6*fs));

% Compute corresponding times & ranges
detected_times = t_mf(locs);
detected_ranges = detected_times * c / 2;

% Identify target as largest peak
[~, max_idx] = max(pks);
target_time = detected_times(max_idx);
target_range = detected_ranges(max_idx);

fprintf('Detected Target Range = %.2f km\n', target_range/1000);

%% ====================== PLOTS ======================
figure('Position',[100 100 900 900]);

% --- Transmitted Rectangular Pulse ---
subplot(3,1,1);
plot(t_tx*1e6, tx_rect, 'LineWidth', 1.5);
title('Transmitted Rectangular Pulse (Long Time Axis)');
xlabel('Time (\mus)');
ylabel('Amplitude');
grid on;

% --- Received Signal ---
subplot(3,1,2);
plot(t_rx*1e6, rx_signal, 'LineWidth', 1);
title('Received Signal (Target + Multiple Clutter + Noise)');
xlabel('Time (\mus)');
ylabel('Amplitude');
grid on;

% --- Matched Filter Output ---
subplot(3,1,3);
plot(t_mf*1e6, abs(mf_output), 'b', 'LineWidth', 1.3);
hold on;
plot(detected_times*1e6, pks, 'ro', 'MarkerFaceColor','r');
title('Matched Filter Output (Pulse Compression)');
xlabel('Time (\mus)');
ylabel('|Amplitude|');
legend('Matched Filter Output','Detected Peaks');
grid on;
% Display all detected ranges
disp('Detected Echo Ranges (km):');
disp(detected_ranges'/1000);
```
