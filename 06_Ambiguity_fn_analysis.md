# Version 1
```matlab
clc; clear; close all;

%% 1. Waveform Parameters (you can change these to test different waveforms)

% Choose waveform type: 'LFM', 'NLFM', 'PhaseCoded', 'Unmodulated'
waveformType = 'LFM';          % Try 'LFM', 'PhaseCoded', or 'Unmodulated'

% Common parameters
fs       = 100e6;              % Sampling frequency (Hz) → high to avoid aliasing
T        = 10e-6;              % Pulse duration (s)
B        = 20e6;               % Bandwidth (Hz) → only used for LFM/NLFM
PRI      = 100e-6;             % Pulse Repetition Interval (s) (not used in AF)

% PhaseCode = 'Frank';         % 'Frank', 'P1', 'P3', 'P4', 'Barker13', etc.

%% 2. Generate Baseband Complex Envelope s(t)

N  = round(T * fs);                    % Number of samples
t  = (0:N-1)' / fs - 0.5/fs;             % Time vector centered at zero
s  = zeros(N,1);                       % Transmit signal (complex envelope)

switch upper(waveformType)
    case 'UNMODULATED'               % Rectangular pulse (constant amplitude)
        s = ones(N,1);

    case 'LFM'                       % Linear Frequency Modulation (Chirp)
        K = B/T;                     % Chirp rate (Hz/s)
        s = exp(1j*pi*K*t.^2);       % Up-chirp

    case 'NLFM'                      % Example of simple Nonlinear FM (S-shaped)
        % Simple tanh-based NLFM for better sidelobes
        alpha = 5;
        instantaneous_freq = B/2 * tanh(alpha * (2*t/T - 1));
        phase = 2*pi * cumtrapz(t, instantaneous_freq);
        s = exp(1j * phase);

    case 'PHASECODED'
        % Example: 13-element Barker code
        if strcmpi(PhaseCode,'Barker13')
            barker13 = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1];
            chips    = repmat(barker13, ceil(N/length(barker13)), 1);
            chips    = chips(1:N);
            s        = chips(:);
        elseif contains(upper(PhaseCode),{'FRANK','P1','P2','P3','P4'})
            M = 8;                  % Change M for different code length (M² samples)
            code = frank_code(M);    % You can implement other codes similarly
            chips = repmat(code, ceil(N/numel(code)), 1);
            s = chips(1:N);
        else
            error('Unknown phase code');
        end

    otherwise
        error('Unknown waveform type');
end

% Normalize energy to 1
s = s / sqrt(s'*s);

%% 3. Compute Ambiguity Function |χ(τ,ν)|

% Delay axis (range)
maxDelaySamples = N-1;
tau = (-maxDelaySamples:maxDelaySamples)' / fs;   % delay in seconds

% Doppler axis – choose reasonable range
fd_max = 100e3;                       % Max Doppler of interest (±100 kHz)
N_doppler = 401;                      % Number of Doppler bins (odd → includes zero)
fd = linspace(-fd_max, fd_max, N_doppler); % Doppler frequency vector (Hz)

AF = zeros(numel(tau), numel(fd));

% Zero-padding the signal for fine delay resolution
s_padded = [zeros(maxDelaySamples,1); s; zeros(maxDelaySamples,1)];

for k = 1:numel(fd)
    % Doppler shift in time domain
    doppler_phase = exp(1j*2*pi*fd(k)*(t - t(1)));
    s_dop = s .* doppler_phase;

    % Correlation (matched filter output) for every delay
    corr = conv(s_padded, conj(flipud(s_dop)));
    % Extract the central part corresponding to our tau vector
    center = length(s_padded);
    AF(:,k) = abs(corr(center-maxDelaySamples:center+maxDelaySamples));
end

% Normalize so that χ(0,0) = 1
AF = AF / max(AF(:));

%% 4. Volume check (should be approximately 1 for any waveform)
volume = sum(abs(AF(:)).^2) * (tau(2)-tau(1)) * (fd(2)-fd(1));
fprintf('Ambiguity function volume = %.6f (should be ≈1)\n', volume);

%% 5. Plotting Results

figure('Color','w','Position',[100 100 1200 900]);

% 5.1 2D Contour plot (classic ambiguity diagram)
subplot(1,2,1);
contourf(fd/1e3, tau*3e8/2/1e3, 20*log10(AF+1e-12), 50, 'LineColor','none');
xlabel('Doppler frequency f_d (kHz)');
ylabel('Range (km)');
title('Ambiguity Function |χ(τ,ν)| (dB)');
cb = colorbar; ylabel(cb,'Magnitude (dB)');
caxis([-60 0]);
colormap jet;
grid on;

% Mark the origin
hold on; plot(0,0,'wp','MarkerSize',10,'MarkerFaceColor','w'); hold off;

% 5.2 3D Surface
subplot(1,2,2);
surf(fd/1e3, tau*3e8/2/1e3, 20*log10(AF+1e-12), 'EdgeColor','none');
xlabel('Doppler (kHz)'); ylabel('Range (km)'); zlabel('dB');
title('3D Ambiguity Function');
view(3); grid on; colormap jet; colorbar; caxis([-60 0]);


%% Optional: Save figures
% print('-dpng', '-r300', sprintf('AmbiguityFunction_%s.png', waveformType));

%% Helper function for Frank code (add at the end of the file)
function code = frank_code(M)
    % Generates Frank polyphase code of length M²
    k = 0:M-1;
    code = zeros(M);
    for m = 1:M
        code(m,:) = exp(1j * 2*pi/M * (m-1) * k);
    end
    code = code(:);
    code = code / abs(code);  % unit magnitude
end
```

# Version 2
```matlab
%======================================================================
% Program Title: Ambiguity Function Analysis for Single Target Radar Detection
% Aim          : To simulate a radar waveform and compute its ambiguity
%                function to analyze range and Doppler characteristics
% Software     : MATLAB R2023a or later recommended
% Author       : (Your Name/Registration Number)
% Date         : November 2025
%======================================================================

clc; clear; close all;

%% 1. Waveform Parameters (you can change these to test different waveforms)

% Choose waveform type: 'LFM', 'NLFM', 'PhaseCoded', 'Unmodulated'
waveformType = 'LFM';          % Try 'LFM', 'PhaseCoded', or 'Unmodulated'

% Common parameters
fs       = 100e6;              % Sampling frequency (Hz) → high to avoid aliasing
T        = 10e-6;              % Pulse duration (s)
B        = 20e6;               % Bandwidth (Hz) → only used for LFM/NLFM
PRI      = 100e-6;             % Pulse Repetition Interval (s) (not used in AF)

% PhaseCode = 'Frank';         % 'Frank', 'P1', 'P3', 'P4', 'Barker13', etc.

%% 2. Generate Baseband Complex Envelope s(t)

N  = round(T * fs);                    % Number of samples
t  = (0:N-1)' / fs - 0.5/fs;             % Time vector centered at zero
s  = zeros(N,1);                       % Transmit signal (complex envelope)

switch upper(waveformType)
    case 'UNMODULATED'               % Rectangular pulse (constant amplitude)
        s = ones(N,1);

    case 'LFM'                       % Linear Frequency Modulation (Chirp)
        K = B/T;                     % Chirp rate (Hz/s)
        s = exp(1j*pi*K*t.^2);       % Up-chirp

    case 'NLFM'                      % Example of simple Nonlinear FM (S-shaped)
        % Simple tanh-based NLFM for better sidelobes
        alpha = 5;
        instantaneous_freq = B/2 * tanh(alpha * (2*t/T - 1));
        phase = 2*pi * cumtrapz(t, instantaneous_freq);
        s = exp(1j * phase);

    case 'PHASECODED'
        % Example: 13-element Barker code
        if strcmpi(PhaseCode,'Barker13')
            barker13 = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1];
            chips    = repmat(barker13, ceil(N/length(barker13)), 1);
            chips    = chips(1:N);
            s        = chips(:);
        elseif contains(upper(PhaseCode),{'FRANK','P1','P2','P3','P4'})
            M = 8;                  % Change M for different code length (M² samples)
            code = frank_code(M);    % You can implement other codes similarly
            chips = repmat(code, ceil(N/numel(code)), 1);
            s = chips(1:N);
        else
            error('Unknown phase code');
        end

    otherwise
        error('Unknown waveform type');
end

% Normalize energy to 1
s = s / sqrt(s'*s);

%% 3. Compute Ambiguity Function |χ(τ,ν)|

% Delay axis (range)
maxDelaySamples = N-1;
tau = (-maxDelaySamples:maxDelaySamples)' / fs;   % delay in seconds

% Doppler axis – choose reasonable range
fd_max = 100e3;                       % Max Doppler of interest (±100 kHz)
N_doppler = 401;                      % Number of Doppler bins (odd → includes zero)
fd = linspace(-fd_max, fd_max, N_doppler); % Doppler frequency vector (Hz)

AF = zeros(numel(tau), numel(fd));

% Zero-padding the signal for fine delay resolution
s_padded = [zeros(maxDelaySamples,1); s; zeros(maxDelaySamples,1)];

for k = 1:numel(fd)
    % Doppler shift in time domain
    doppler_phase = exp(1j*2*pi*fd(k)*(t - t(1)));
    s_dop = s .* doppler_phase;

    % Correlation (matched filter output) for every delay
    corr = conv(s_padded, conj(flipud(s_dop)));
    % Extract the central part corresponding to our tau vector
    center = length(s_padded);
    AF(:,k) = abs(corr(center-maxDelaySamples:center+maxDelaySamples));
end

% Normalize so that χ(0,0) = 1
AF = AF / max(AF(:));

%% 4. Volume check (should be approximately 1 for any waveform)
volume = sum(abs(AF(:)).^2) * (tau(2)-tau(1)) * (fd(2)-fd(1));
fprintf('Ambiguity function volume = %.6f (should be ≈1)\n', volume);

%% 5. Plotting Results

figure('Color','w','Position',[100 100 1200 900]);

% 5.1 2D Contour plot (classic ambiguity diagram)
subplot(2,2,1);
contourf(fd/1e3, tau*3e8/2/1e3, 20*log10(AF+1e-12), 50, 'LineColor','none');
xlabel('Doppler frequency f_d (kHz)');
ylabel('Range (km)');
title('Ambiguity Function |χ(τ,ν)| (dB)');
cb = colorbar; ylabel(cb,'Magnitude (dB)');
caxis([-60 0]);
colormap jet;
grid on;

% Mark the origin
hold on; plot(0,0,'wp','MarkerSize',10,'MarkerFaceColor','w'); hold off;

% 5.2 3D Surface
subplot(2,2,2);
surf(fd/1e3, tau*3e8/2/1e3, 20*log10(AF+1e-12), 'EdgeColor','none');
xlabel('Doppler (kHz)'); ylabel('Range (km)'); zlabel('dB');
title('3D Ambiguity Function');
view(3); grid on; colormap jet; colorbar; caxis([-60 0]);

% 5.3 Zero-Doppler cut (autocorrelation → range sidelobes)
subplot(2,2,3);
zeroDopIdx = find(fd==0);
plot(tau*3e8/2/1e3, 20*log10(abs(AF(:,zeroDopIdx))+1e-12), 'LineWidth',1.5);
xlabel('Range (km)'); ylabel('Magnitude (dB)');
title('Zero-Doppler Cut (Range Profile)');
grid on; xlim([-2 2]); ylim([-80 5]);

% 5.4 Zero-delay cut (Doppler tolerance)
subplot(2,2,4);
zeroDelayIdx = find(abs(tau)<1/fs,1);  % closest to τ=0
plot(fd/1e3, 20*log10(abs(AF(zeroDelayIdx,:))+1e-12), 'LineWidth',1.5);
xlabel('Doppler frequency (kHz)'); ylabel('Magnitude (dB)');
title('Zero-Delay Cut (Doppler Profile)');
grid on; xlim([-fd_max/1e3 fd_max/1e3]); ylim([-80 5]);

sgtitle(sprintf('Ambiguity Function – %s Waveform\nT=%.1f µs, B=%.1f MHz', ...
    waveformType, T*1e6, B/1e6), 'FontSize',14,'FontWeight','bold');

%% Optional: Save figures
% print('-dpng', '-r300', sprintf('AmbiguityFunction_%s.png', waveformType));

%% Helper function for Frank code (add at the end of the file)
function code = frank_code(M)
    % Generates Frank polyphase code of length M²
    k = 0:M-1;
    code = zeros(M);
    for m = 1:M
        code(m,:) = exp(1j * 2*pi/M * (m-1) * k);
    end
    code = code(:);
    code = code / abs(code);  % unit magnitude
end
```

# Version 3 : final
```matlab
% Matched Filter Output / Cross-Ambiguity Function for Multiple Targets
clear; clc; close all;

%% ============== USER / DEFAULT SIGNAL PARAMETERS ==============
% Basic transmitted/received complex samples (row vector)
% Using a high-resolution LFM (Chirp) for good separation
L_basic = 16;
u_basic = exp(1j * pi * (0:L_basic-1).^2 / L_basic);

% Oversampling ratio (samples per tb)
r = 4;

% Number of Doppler bins (K)
K = 512;

% Max integer delay (in sample-units after oversampling). delays = -N..N
N = 120;
% =====================================================================

%% ============== MULTI-TARGET DEFINITION (change here) ==============
% Define targets as rows: [Amplitude, Delay (in tb units), Doppler (cycles/tb)]
% NOTE: Delay and Doppler are relative to the grid resolution.
targets = [
    1.0, 0.0, 0.0;     % Target 1: Strong, reference (near origin)
    0.8, 10.0, 0.4;    % Target 2: Moderate, delayed, positive Doppler
    0.5, -5.0, -0.2;   % Target 3: Weak, early, negative Doppler
    0.6, 25.0, 0.1;    % Target 4: Distant, small positive Doppler
];

% Additive White Gaussian Noise power (relative to signal peak)
SNR_dB = 20; % 20 dB SNR
% =====================================================================

% ensure row vector
u_basic = u_basic(:).';

% Simple repeat upsampling (each basic sample repeated r times)
u = kron(u_basic, ones(1, r));
L = length(u);
fs = r; % sampling frequency in samples per tb

% ------------------- Matched Filter / Cross-Ambiguity Setup -------------------
% delays (integer sample shifts)
delays = -N : N;
numDelays = numel(delays);

% Doppler axis (cycles per tb)
fd = ((-K/2 : K/2-1) / K) * fs;% row vector length K

% physical axes
tau = delays / fs;% delay in units of tb
dopp = fd;% Doppler in cycles per tb

% ------------------- 1. Construct Received Signal r(t) -------------------
T_samp = 1/fs; % sample time
t_idx = 0 : L-1; % time indices
t_vec = t_idx * T_samp; % time vector (in tb units)

r_total = zeros(1, L);

% Calculate the maximum signal power for noise scaling
max_signal_power = max(abs(u).^2);
noise_power = max_signal_power / (10^(SNR_dB/10));
std_dev = sqrt(noise_power / 2); % Divide by 2 for real and imag parts

% Loop through all targets
for i = 1:size(targets, 1)
    A_k = targets(i, 1);
    tau_k = targets(i, 2);
    fd_k = targets(i, 3);
   
    % The integer delay (in samples) for this target
    d_samp = round(tau_k * fs);
   
    % The index shift in the array
    shift = -d_samp;
   
    % Define the target's complex return (Doppler shift)
    doppler_term = exp(1j * 2 * pi * fd_k * t_vec);
    u_doppled = u .* doppler_term;
   
    % Shift the doppled signal to create the delayed return
    r_k = zeros(1, L);
   
    if shift >= 0 % Target delayed (right shift, positive tau)
        % u(1:L-shift) goes into r_k(1+shift:L)
        r_k(1+shift : L) = u_doppled(1 : L-shift);
    else % Target advanced (left shift, negative tau)
        % u(1-shift:L) goes into r_k(1 : L+shift)
        r_k(1 : L+shift) = u_doppled(1-shift : L);
    end
   
    % Sum all returns
    r_total = r_total + A_k * r_k;
end

% Add AWGN
noise = std_dev * (randn(1, L) + 1j * randn(1, L));
r_total = r_total + noise;


% ------------------- 2. Compute Cross-Ambiguity Function -------------------
A_cross = zeros(K, numDelays);% complex cross-ambiguity (Doppler x delay)

% Matched filter loop: Cross-correlation of r(t) with u(t)
for col = 1 : numDelays
    d = delays(col); % Integer delay index (tau)
   
    % A_cross(d, fd) = sum_t r(t) * conj(u(t-d)) * exp(-j 2 pi fd t)
   
    if d >= 0
        % r(t) is r(1+d:L), u(t-d) is u(1:L-d)
        len = L - d;
        if len > 0
            prod = r_total(1 + d : d + len) .* conj(u(1 : len)); % r(t) * conj(u(t-d))
        else
            prod = zeros(1,0);
        end
    else % d < 0
        dpos = -d;
        len = L - dpos;
        if len > 0
            % r(t) is r(1:L-dpos), u(t-d) is u(1+dpos:L)
            prod = r_total(1 : len) .* conj(u(1 + dpos : dpos + len)); % r(t) * conj(u(t-d))
        else
            prod = zeros(1,0);
        end
    end

    if isempty(prod)
        A_cross(:, col) = zeros(K,1);
    else
        prod = prod(:);
        S = fft(prod, K);
        S = fftshift(S);
        A_cross(:, col) = S;
    end
end

% Magnitude and normalization
A_mag = abs(A_cross);
maxA = max(A_mag(:));
if maxA == 0
    warning('Cross-Ambiguity magnitude is zero everywhere.');
    A_norm = A_mag;
else
    A_norm = A_mag / maxA;% normalized to peak 1
end

% ------------------- 3. Plotting Results -------------------

%% ---- PLOT 1: 2D contour (dB) to clearly show target peaks ----
figure(1); clf;
epsdb = 1e-12;
AdB = 20*log10(A_norm + epsdb);
contourf(tau, dopp, AdB, 40, 'LineColor','none');
hold on;

% Plot markers for true target positions
for i = 1:size(targets, 1)
    plot(targets(i, 2), targets(i, 3), 'wx', 'MarkerSize', 10, 'LineWidth', 2);
end

xlabel('\tau (tb units)','FontSize',12);
ylabel('f_d (cycles / tb)','FontSize',12);
title('Figure 1: Cross-Ambiguity Function |A_{r,u}| (dB)','FontSize',14);
colorbar;
caxis([-60 0]); % Adjust dynamic range to highlight targets
grid on;

%% ---- PLOT 2: 3D surface (magnitude) ----
figure(2); clf;
h = surf(tau, dopp, A_norm, 'EdgeColor','none');
axis tight;
view(40,30);
xlabel('\tau (tb units)','FontSize',12);
ylabel('f_d (cycles / tb)','FontSize',12);
zlabel('|A_{r,u}(\tau,f_d)| (normalized)','FontSize',12);
title('Figure 2: Cross-Ambiguity Function (magnitude)','FontSize',14);
colormap(jet);
colorbar;
set(h,'FaceLighting','phong','AmbientStrength',0.3);

%% ---- PLOT 3: Important Cuts (Autocorrelation and Spectrum) ----
% Find the indices corresponding to the main target (Target 1: tau=0, fd=0)
[~, idx_fd0] = min(abs(dopp));
[~, idx_tau0] = min(abs(tau));

figure(3); clf;

% Subplot 1: Delay Cut (Autocorrelation)
subplot(2, 1, 1);
plot(tau, A_norm(idx_fd0, :), 'LineWidth', 2);
title('Figure 3a: Delay Cut |A_{r,u}(\tau, f_d=0)| - Highlights Target 2 and 4');
xlabel('\tau (tb units)');
ylabel('|A| (normalized)');
grid on;

% Subplot 2: Doppler Cut (Spectrum)
subplot(2, 1, 2);
plot(dopp, A_norm(:, idx_tau0), 'LineWidth', 2);
title('Figure 3b: Doppler Cut |A_{r,u}(\tau=0, f_d)| - Highlights Target 1');
xlabel('f_d (cycles / tb)');
ylabel('|A| (normalized)');
grid on;

disp('Multi-target Cross-Ambiguity computation finished.');
disp('Figure 1 shows the 2D contour with white Xs marking the true target locations.');
disp('Figure 3 cuts are taken through the strongest target at (tau=0, fd=0).');
```
