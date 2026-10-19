% PPG HEART RATE - 53 SUBJECTS

load('C:\Users\barad\OneDrive\Documents\MATLAB\dsp\bidmc_data.mat')
numSubjects = length(data);
Results = zeros(numSubjects, 7);
ecg_prom_factor = 2.0;

%loop to run on each subject
for s = 1:numSubjects
    fprintf('Processing Subject %02d ... ', s)
    try
        
        %Pre-prcoessing        
        % Extract raw signal vectors and sampling properties
        ppg = data(s).ppg.v;
        Fs  = data(s).ppg.fs;

        % Center around zero to strip out massive DC voltage offsets
        ppg_dc = ppg - mean(ppg);

        % Apply 4th-order Butterworth bandpass (0.5 to 5 Hz / 30 to 300 BPM)
        lowCut  = 0.5;
        highCut = 5;
        [b, a] = butter(4, [lowCut highCut] / (Fs/2), 'bandpass');
        ppg_bw = filtfilt(b, a, ppg_dc);

        % Smooth local high-frequency digital ripples
        ppg_ma = movmean(ppg_bw, 5);

        % Restore peak shapes while polishing out remaining noise
        ppg_final = sgolayfilt(ppg_ma, 3, 11);
        
        % Find heartbeat peaks with a 0.4s refractory lockout period
        [pks, locs] = findpeaks(ppg_final, ...
            'MinPeakDistance', round(0.4 * Fs), ...
            'MinPeakProminence', 0.05);

        % Translate sample index spaces into absolute physical time
        peak_time = locs / Fs;

        % Calculate peak-to-peak durations (RR Intervals)
        RR_interval = diff(peak_time);

        % Filter intervals to strictly human ranges (0.3s to 2.0s)
        RR_valid = RR_interval(RR_interval > 0.3 & RR_interval < 2.0);

        % Convert validated durations directly into Beats Per Minute
        HR = 60 ./ RR_valid;
        average_HR = mean(HR);

        %using ecg for reference
        
        % Extract reference cardiovascular metrics
        ecg    = data(s).ekg.v;
        Fs_ecg = data(s).ekg.fs;
        ecg_dc = ecg - mean(ecg);

        % Compute an adaptive peak threshold unique to this subject's variance
        ecg_prominence = ecg_prom_factor * std(ecg_dc);

        % Find reference electrical R-wave spikes
        [ecg_peaks, ecg_locs] = findpeaks(ecg_dc, ...
            'MinPeakDistance', round(0.4 * Fs_ecg), ...
            'MinPeakProminence', ecg_prominence);

        % Convert sample indexes to absolute timelines
        ecg_peak_time = ecg_locs / Fs_ecg;
        ecg_RR = diff(ecg_peak_time);

        % Validate reference intervals against physiological bounds
        ecg_RR_valid = ecg_RR(ecg_RR > 0.3 & ecg_RR < 2.0);

        % Calculate the definitive reference Heart Rate
        ecg_HR = 60 ./ ecg_RR_valid;
        reference_HR = mean(ecg_HR);

        %statistical analysis
        
        % Measure divergence between optical PPG and electrical ECG systems
        HR_error         = abs(average_HR - reference_HR);
        HR_percent_error = (HR_error / reference_HR) * 100;

        % Commit telemetry array records to memory
        Results(s, :) = [ ...
            s, ...
            average_HR, ...
            reference_HR, ...
            HR_error, ...
            HR_percent_error, ...
            length(locs), ...
            length(RR_valid)];

        if isnan(reference_HR)
            fprintf('SUCCESS (ECG still unresolved)\n')
        else
            fprintf('SUCCESS\n')
        end

    catch ME
        % Handle unexpected signal drops or missing data vectors gracefully
        fprintf('FAILED\n')
        fprintf('   %s\n', ME.message)

        Results(s, :) = [s, NaN, NaN, NaN, NaN, NaN, NaN];
    end
end

%display findings
% Construct an structured analysis table for workspace investigation
ResultsTable = array2table(Results, ...
    'VariableNames', { ...
    'Subject', ...
    'PPG_HR', ...
    'ECG_HR', ...
    'HR_Error', ...
    'HR_Percent_Error', ...
    'Detected_Beats', ...
    'Valid_RR_Intervals'});

fprintf('RESULTS\n')
disp(ResultsTable)

% Isolate completely clean tracking runs via error validity checks
valid = ~isnan(Results(:, 4));
fprintf('Total subjects        = %d\n', numSubjects)
fprintf('Valid subjects        = %d\n', sum(valid))

if sum(valid) < numSubjects
    stillBad = Results(~valid, 1);
    fprintf('Still unresolved      = Subject(s) %s\n', mat2str(stillBad'))
end

fprintf('\nMean PPG HR           = %.2f BPM\n', mean(Results(:, 2), 'omitnan'))
fprintf('Mean ECG HR           = %.2f BPM\n', mean(Results(:, 3), 'omitnan'))
fprintf('Mean Abs Error        = %.2f BPM\n', mean(Results(:, 4), 'omitnan'))
fprintf('Mean Percent Error    = %.2f%%\n',  mean(Results(:, 5), 'omitnan'))
