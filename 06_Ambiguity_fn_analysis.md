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
