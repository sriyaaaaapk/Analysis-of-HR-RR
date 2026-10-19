%for visualization only

clc
clear
close all
%HR PIPELINE VISUALIZATION for SUBJECT 1
load('C:\Users\barad\OneDrive\Documents\MATLAB\dsp\bidmc_data.mat')
subject = data(1);

ppg = subject.ppg.v;
Fs  = subject.ppg.fs;

t = (0:length(ppg)-1) / Fs;
fprintf('Sampling Frequency = %.2f Hz\n', Fs)
fprintf('Signal Duration    = %.2f seconds\n', length(ppg)/Fs)


%Use only first 20 Seconds
idx = t <= 20;

% STEP 1 - RAW PPG
ppg_raw = ppg;

%STEP 2 - DC REMOVAL
ppg_dc = ppg_raw - mean(ppg_raw);

% Step 3 - BUTTERWORTH : set as per human freq of HR
lowCut  = 0.5;
highCut = 5;
[b, a] = butter(4, ...
    [lowCut highCut] / (Fs/2), ...
    'bandpass');

ppg_bw = filtfilt(b, a, ppg_dc);

%STEP 4 - MOVING AVERAGE
windowSize = 5;
ppg_ma = movmean(ppg_bw, windowSize);

%STEP 5 - SAVITZKY-GOLAY
polyOrder = 3;
frameSize = 11;

ppg_final = sgolayfilt(ppg_ma, polyOrder, frameSize);

%STEP 6 - PEAK DETECTION
[pks, locs] = findpeaks(ppg_final, ...
    'MinPeakDistance', round(0.4 * Fs), ...
    'MinPeakProminence', 0.05);
peak_idx = locs(t(locs) <= 20);

%figures
figure('Name', 'HR Pipeline - Filtering', ...
       'Position', [100 50 1300 900])
subplot(5,1,1)
plot(t(idx), ppg_raw(idx), 'LineWidth', 1)
xlabel('Time (seconds)')
ylabel('Amplitude')
title('1. Raw PPG Signal')
grid on

subplot(5,1,2)
plot(t(idx), ppg_dc(idx), 'LineWidth', 1)
xlabel('Time (seconds)')
ylabel('Amplitude')
title('2. DC Offset Removed')
grid on

subplot(5,1,3)
plot(t(idx), ppg_bw(idx), 'LineWidth', 1)
xlabel('Time (seconds)')
ylabel('Amplitude')
title('3. 4th-Order Butterworth Bandpass (0.5–5 Hz)')
grid on

subplot(5,1,4)
plot(t(idx), ppg_ma(idx), 'LineWidth', 1)
xlabel('Time (seconds)')
ylabel('Amplitude')
title('4. Moving Average Filter (Window = 5)')
grid on

subplot(5,1,5)
plot(t(idx), ppg_final(idx), 'LineWidth', 1)
xlabel('Time (seconds)')
ylabel('Amplitude')
title('5. Savitzky-Golay Filter (Order = 3, Frame = 11)')
grid on
sgtitle('PPG Heart Rate Processing Pipeline')

%peak detection
figure('Name', 'HR Pipeline - Peak Detection', ...
       'Position', [150 100 1300 500])
plot(t(idx), ppg_final(idx), 'LineWidth', 1)
hold on
plot(t(peak_idx), ...
     ppg_final(peak_idx), ...
     'ro', ...
     'MarkerFaceColor', 'r')
xlabel('Time (seconds)')
ylabel('Amplitude')
title('Final PPG Signal with Detected Heartbeat Peaks')
legend('Processed PPG', 'Detected Peaks')
grid on
