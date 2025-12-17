# Expt-10 : Adaptive Radar Target Detection Using CFAR

## Version-1
```matlab
%==========================================================================
% Title   : Adaptive Radar Target Detection using CA-CFAR
% Aim     : To simulate adaptive radar detection using Cell-Averaging CFAR
% Software: MATLAB
%==========================================================================

clc;
clear;
close all;

%% 1. Radar and CFAR Parameters

numCells     = 1000;        % Total number of range cells
numTargets   = 3;           % Number of targets
targetPos    = [200 500 750]; % Target locations (range bins)
targetAmp    = [15 20 18];  % Target amplitudes

Pfa          = 1e-4;        % Desired probability of false alarm
numTrain     = 20;          % Training cells on each side
numGuard     = 4;           % Guard cells on each side

%% 2. Generate Noise (Complex Gaussian)

noise = (randn(1,numCells) + 1j*randn(1,numCells))/sqrt(2);
signal = abs(noise).^2;     % Power of received signal

%% 3. Insert Targets

for k = 1:numTargets
    signal(targetPos(k)) = signal(targetPos(k)) + targetAmp(k)^2;
end

%% 4. CA-CFAR Processing

N = 2 * numTrain;   % Total training cells
alpha = N * (Pfa^(-1/N) - 1);   % CFAR scaling factor

threshold = zeros(1,numCells);
detection = zeros(1,numCells);

for i = numTrain+numGuard+1 : numCells-(numTrain+numGuard)
    
    % Training cells
    leadingCells = signal(i-numGuard-numTrain : i-numGuard-1);
    laggingCells = signal(i+numGuard+1 : i+numGuard+numTrain);
    
    noiseEstimate = mean([leadingCells laggingCells]);
    
    threshold(i) = alpha * noiseEstimate;
    
    % Detection decision
    if signal(i) > threshold(i)
        detection(i) = signal(i);
    end
end

%% 5. Plot Results

figure;
plot(signal,'b','LineWidth',1.2); hold on;
plot(threshold,'r--','LineWidth',1.5);
stem(detection,'g','filled');

xlabel('Range Cell Index');
ylabel('Power');
title('Adaptive Radar Target Detection using CA-CFAR');
legend('Received Signal','CFAR Threshold','Detected Targets');
grid on;
```
