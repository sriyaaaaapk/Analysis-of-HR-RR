%% hardware.m  --  MAX30102 (via Arduino) -> MATLAB: heart rate + respiratory rate from a captured PPG segment
% NOT real-time: it records for hw.duration_s seconds, analyses it, then draws ONE figure and prints HR and RR.
%
% ARDUINO SIDE: only requireement = ONE sample per line, IR value in a known column. Examples of valid lines:
%     52314                      (IR only)                         -> hw.col = 1
%     183422,52314               (millis(), IR)                    -> hw.time_col = 1, hw.col = 2
%     183422,52314,48101         (millis(), IR, RED)               -> hw.time_col = 1, hw.col = 2
%   Output rate = sampleRate / sampleAverage (SparkFun library). E.g. sampleRate 400 & sampleAverage 4 -> 100 Hz. Aim for >= 50 Hz.
%   Close the Arduino Serial Monitor / Plotter first: only one program can own the COM port.
%
% PIPELINE (reuses the project's own functions; put the lin_rqipca folder on the MATLAB path):
%   HR : DC removal -> band-pass 0.5-5 Hz (Butterworth) -> Savitzky-Golay -> peaks -> beat intervals -> HR = 60 / interval
%   RR : resample to 300 Hz -> prep_signal -> 30 s windows -> PPG components PAM, PFM, PBW (ppg_fiducials) -> 4 Hz, z-score
%        -> vmd_refine -> select components ('nb' = classifier trained on BIDMC, 'rqi' = Lin rule, 'all') -> pca_fuse -> RR = 60*f_peak
% NEEDS: config, prep_signal, ppg_fiducials, vmd_refine, pca_fuse, rr_from_signal, compute_rqi, dom_freq, spec1 (lin_rqipca),
%        Signal Processing Toolbox. For hw.select = 'nb' also results_rmc_cache.mat (made by rmc_classifier_selection.m).
%
% HONEST LIMITS: RR from a finger PPG alone is the weakest setting (on BIDMC the PPG components PAM/PFM/PBW were correct alone in only
%   22-60% of windows), and the selector was trained on BIDMC ICU data, not on this sensor. Treat RR as an estimate; check it against a
%   manual 30-60 s breath count (hw.manual_rr) while you sit still, breathe normally and keep the finger steady.

clear; clc; close all;
%% ---------------- SETTINGS ----------------
hw.source     = 'serial';      % 'serial' (live capture) | 'file' (re-analyse a saved capture) | 'demo' (synthetic test signal)
hw.port       = 'COM5';        % your Arduino port; the list of available ports is printed at start
hw.baud       = 115200;
hw.duration_s = 90;            % capture length: >= 60 s recommended (RR uses 30 s windows)
hw.col        = 1;             % column holding the IR value in each line
hw.time_col   = 0;             % column holding millis() (0 = none: sampling rate is estimated from the capture time)
hw.fs         = [];            % [] = automatic; or force e.g. 100
hw.file       = 'capture.csv'; % used when hw.source = 'file'; live captures are saved to a time-stamped CSV
hw.select     = 'nb';          % RR component selection: 'nb' | 'rqi' | 'all'
hw.cache      = 'results_rmc_cache.mat';
hw.invert     = 'auto';        % PPG polarity: 'auto' | true | false (reflective finger sensors are often upside-down)
hw.plot_start_s = 10; hw.plot_len_s = 12;   % PPG segment shown in panels 1-2
hw.rr_window  = 'last';        % which 30 s window is drawn in panels 4-5: 'last' or an index
hw.manual_rr  = NaN;           % optional: breaths/min you counted by hand (shown in the figure)
hw.manual_hr  = NaN;           % optional: reference heart rate (e.g. from a watch)

%% ---------------- 1. GET RAW PPG ----------------
switch hw.source
    case 'serial'
        fprintf('Available ports: %s\n', strjoin(string(serialportlist("available")), ', '));
        s = serialport(hw.port, hw.baud); configureTerminator(s, 'LF'); flush(s); pause(0.5); flush(s);
        maxN = hw.duration_s*600; X = nan(maxN, 3); tt = nan(maxN, 1); k = 0; t0 = tic;
        fprintf('Recording %d s ... keep the finger still, breathe normally.\n', hw.duration_s);
        try
            while toc(t0) < hw.duration_s
                if s.NumBytesAvailable > 0
                    v = str2double(regexp(char(readline(s)), '[-+]?\d+\.?\d*', 'match'));
                    if numel(v) >= max(hw.col, hw.time_col) && k < maxN
                        k = k + 1; v = v(1:min(3, numel(v))); X(k, 1:numel(v)) = v; tt(k) = toc(t0);
                    end
                else, pause(0.001);
                end
            end
        catch ME, clear s; rethrow(ME);
        end
        clear s; X = X(11:k, :); tt = tt(11:k);                      % drop first samples (possible partial lines)
        fname = ['hardware_capture_' datestr(now, 'yyyymmdd_HHMMSS') '.csv']; writematrix(X(:, any(isfinite(X), 1)), fname);
        fprintf('Saved %d samples to %s\n', size(X,1), fname);
    case 'file'
        X = readmatrix(hw.file); tt = [];
    case 'demo'                                                      % synthetic PPG: HR 72 BPM, RR 15 br/min, upside-down like a reflective sensor
        fsd = 100; td = (0:1/fsd:hw.duration_s)'; rng(1);
        ph = 2*pi*cumsum(1.2*(1 + 0.05*sin(2*pi*0.25*td)))/fsd; pp = mod(ph/(2*pi), 1);
        pulse = (pp < 0.15).*(pp/0.15) + (pp >= 0.15).*exp(-(pp-0.15)/0.25);
        sg = pulse.*(1 + 0.15*sin(2*pi*0.25*td)) + 0.1*sin(2*pi*0.25*td + 1);
        X = [td*1000, 100000 - 2000*sg + 20*randn(size(td))]; hw.time_col = 1; hw.col = 2; tt = [];
end
raw = X(:, hw.col);
if hw.time_col > 0                                                   % timestamps -> uniform grid
    tms = X(:, hw.time_col)/1000; fs = round(1/median(diff(tms)));
    if ~isempty(hw.fs), fs = hw.fs; end
    tu = (tms(1):1/fs:tms(end))'; raw = interp1(tms, raw, tu, 'linear');
else
    if ~isempty(hw.fs), fs = hw.fs; elseif ~isempty(tt), fs = round((numel(tt)-1)/(tt(end)-tt(1)));
    else, error('No timestamps in the file: set hw.fs.'); end
end
raw = fillmissing(raw(:), 'nearest'); rawU = raw; T = numel(raw)/fs;
fprintf('%d samples, fs = %d Hz, duration %.1f s, IR range %.0f-%.0f\n', numel(raw), fs, T, min(raw), max(raw));
assert(fs >= 20, 'Sampling rate %d Hz is too low (need >= 20, ideally >= 50). Change the Arduino sample settings.', fs);
assert(T >= 35, 'Need at least ~35 s of data for one 30 s RR window (got %.0f s).', T);

%% ---------------- 2. COMMON PRE-PROCESSING ----------------
cfg = config();                                                      % Lin parameters (fs = 300 Hz, VMD 5 modes, ...)
xr = raw - mean(raw);                                                % remove DC BEFORE resampling (avoids edge transients)
[pq, qq] = rat(cfg.fs/fs, 1e-6); x = resample(xr, pq, qq); t = (0:numel(x)-1)'/cfg.fs;
[z, p_, k_] = butter(2, [0.5 5]/(cfg.fs/2), 'bandpass'); [sos, g] = zp2sos(z, p_, k_);
fl = 2*floor(0.06*cfg.fs) + 1;                                       % ~0.12 s Savitzky-Golay frame (odd)
if islogical(hw.invert) || isnumeric(hw.invert), doflip = logical(hw.invert);
else                                                                 % auto: PPG upstroke is faster than the decay
    dd = diff(sgolayfilt(filtfilt(sos, g, x), 3, fl)); c = round(3*cfg.fs):numel(dd)-round(3*cfg.fs);
    doflip = max(dd(c)) < -min(dd(c));
end
if doflip, x = -x; fprintf('Polarity: signal inverted (auto-detected).\n'); else, fprintf('Polarity: signal used as is.\n'); end

%% ---------------- 3. HEART RATE ----------------
xs = sgolayfilt(filtfilt(sos, g, x), 3, fl);                         % "proper" PPG: band-pass 0.5-5 Hz + Savitzky-Golay
[~, loc] = findpeaks(xs, 'MinPeakDistance', round(0.33*cfg.fs), 'MinPeakProminence', 0.4*std(xs));
tpk = t(loc); ibi = diff(tpk); good = ibi > 0.33 & ibi < 1.5;        % 40-180 BPM
hr_inst = 60 ./ ibi(good); t_hr = tpk([false; good]);
HR = 60 / median(ibi(good));
fprintf('\nHEART RATE: %.1f BPM  (median of %d beat intervals; mean %.1f, SD %.1f)\n', HR, nnz(good), mean(hr_inst), std(hr_inst));

%% ---------------- 4. RESPIRATORY RATE (PPG-only components) ----------------
ppg = prep_signal(x, cfg.fs, cfg);                                   % fs_in == cfg.fs -> no resampling; 10-pt smoothing + 0.05 Hz high-pass
Nw = 30*cfg.fs; nW = floor(numel(ppg)/Nw); L = 30*cfg.fs_rmc; tg = (0:L-1)/cfg.fs_rmc; names = {'PAM','PFM','PBW'};
useNB = strcmp(hw.select, 'nb');
if useNB && ~isfile(hw.cache), warning('%s not found: falling back to the RQI rule.', hw.cache); useNB = false; hw.select = 'rqi'; end
if useNB                                                             % per-component naive Bayes trained on BIDMC (components 4-6 = PAM, PFM, PBW)
    S = load(hw.cache, 'C'); Cb = S.C; okb = Cb.ok;
    for i = 1:3, ytr{i} = double(abs(Cb.rr(okb, i+3) - Cb.ref(okb)) <= 2); Xtr{i} = squeeze(Cb.td(okb, i+3, :)); end %#ok<SAGROW>
end
RRw = nan(nW, 1); win = cell(nW, 1);
for w = 1:nW
    seg = ppg((w-1)*Nw+1 : w*Nw); tw = (0:Nw-1)'/cfg.fs;
    [tr, pk] = ppg_fiducials(seg, cfg.fs, cfg.ppg_peak_win);
    if numel(pk) < cfg.min_beats, fprintf('window %d: too few beats, skipped\n', w); continue; end
    R = [grid_(tw(pk), seg(pk)-seg(tr), tg); grid_(tw(pk(2:end)), diff(tw(pk)), tg); grid_(tw(pk), (seg(pk)+seg(tr))/2, tg)];
    Rr = zeros(3, L); rr_i = zeros(1, 3); rq = zeros(1, 3); td = zeros(3, 6);
    for i = 1:3
        Rr(i,:) = vmd_refine(R(i,:), cfg); rr_i(i) = rr_from_signal(Rr(i,:), cfg.fs_rmc, cfg);
        rq(i) = compute_rqi(Rr(i,:), cfg.fs_rmc, cfg); td(i,:) = td_features_(Rr(i,:), 0.01);
    end
    switch hw.select
        case 'nb',  sc = zeros(1, 3); for i = 1:3, [~, q_] = gnb_(Xtr{i}, ytr{i}, td(i,:)); sc(i) = q_; end
        case 'rqi', sc = rq;
        otherwise,  sc = ones(1, 3);
    end
    sel = sc >= mean(sc) - eps; fused = pca_fuse(Rr(sel,:)); RRw(w) = rr_from_signal(fused, cfg.fs_rmc, cfg);
    win{w} = struct('R', R, 'Rr', Rr, 'rr_i', rr_i, 'rq', rq, 'sc', sc, 'sel', sel, 'fused', fused);
    fprintf('window %d (%2d-%2d s): RR = %5.1f br/min | components %s -> own rates %s | kept: %s\n', w, (w-1)*30, w*30, RRw(w), ...
        strjoin(names, '/'), mat2str(round(rr_i, 1)), strjoin(names(sel), '+'));
end
assert(any(isfinite(RRw)), 'No window could be analysed (too few beats?). Check finger placement / polarity / sampling rate.');
RR = median(RRw(isfinite(RRw)));
fprintf('\nRESPIRATORY RATE: %.1f br/min  (median of %d windows, selection = %s)\n', RR, nnz(isfinite(RRw)), hw.select);
if isfinite(hw.manual_rr), fprintf('  manual count %.1f -> error %+.1f br/min\n', hw.manual_rr, RR - hw.manual_rr); end
if isfinite(hw.manual_hr), fprintf('  reference HR %.1f -> error %+.1f BPM\n', hw.manual_hr, HR - hw.manual_hr); end

%% ---------------- 5. ONE FIGURE (a chosen segment, not live) ----------------
wi = hw.rr_window; if ischar(wi) || isstring(wi), wi = find(isfinite(RRw), 1, 'last'); end
if ~isfinite(RRw(wi)), wi = find(isfinite(RRw), 1, 'last'); end
ps = min(hw.plot_start_s, max(0, T - hw.plot_len_s)); seg_t = [ps, ps + hw.plot_len_s];
f = figure('Position', [60 40 1200 900], 'Color', 'w'); tl = tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile([1 2]); tr_ = (0:numel(rawU)-1)'/fs; m = tr_ >= seg_t(1) & tr_ <= seg_t(2); plot(tr_(m), rawU(m), 'Color', [0.1 0.4 0.7], 'LineWidth', 1.2);
ylabel('IR (ADC counts)'); title('1. Raw PPG from the MAX30102 (selected segment)'); grid on; box off; xlim(seg_t);
nexttile([1 2]); m = t >= seg_t(1) & t <= seg_t(2); plot(t(m), xs(m), 'Color', [0.1 0.4 0.7], 'LineWidth', 1.4); hold on
mp = tpk >= seg_t(1) & tpk <= seg_t(2); plot(tpk(mp), xs(loc(mp)), 'v', 'MarkerFaceColor', [0.9 0.2 0.2], 'MarkerEdgeColor', 'w', 'MarkerSize', 8);
ylabel('filtered PPG (a.u.)'); title('2. Proper PPG: band-pass 0.5-5 Hz + Savitzky-Golay, detected beats marked'); grid on; box off; xlim(seg_t);
nexttile([1 2]); plot(t_hr, hr_inst, '-o', 'Color', [0.95 0.5 0.1], 'MarkerFaceColor', [0.95 0.5 0.1], 'MarkerSize', 4); hold on; yline(HR, '--k', sprintf('median %.1f BPM', HR));
ylabel('heart rate (BPM)'); xlabel('time (s)'); title('3. Instantaneous heart rate over the whole recording'); grid on; box off; xlim([0 T]);
nexttile; d = win{wi}; tz = (0:L-1)/cfg.fs_rmc + (wi-1)*30; plot(tz, d.Rr', 'Color', [0.8 0.8 0.8]); hold on
plot(tz, zscore_(d.fused), 'Color', [0.18 0.55 0.43], 'LineWidth', 2); xlabel('time (s)'); ylabel('a.u.');
title(sprintf('4. Respiratory waveform, window %d (grey = components, green = fused)', wi), 'FontSize', 10); grid on; box off; xlim([tz(1) tz(end)]);
nexttile; [fq, Sp] = spec1(d.fused, cfg.fs_rmc, cfg.nfft); fb = fq*60; Sp = Sp/max(Sp(fq >= cfg.fmin & fq <= cfg.fmax));
patch([3 45 45 3], [0 0 1.2 1.2], [0.9 0.9 0.9], 'EdgeColor', 'none'); hold on; plot(fb, Sp, 'Color', [0.18 0.55 0.43], 'LineWidth', 2);
xline(RRw(wi), '--k', sprintf('%.1f br/min', RRw(wi))); xlim([0 60]); ylim([0 1.2]); xlabel('breaths per minute'); ylabel('normalised power');
title('5. Respiratory spectrum (shaded = 3-45 br/min)', 'FontSize', 10); grid on; box off;
txt = sprintf('HR = %.0f BPM     |     RR = %.1f br/min', HR, RR);
if isfinite(hw.manual_rr), txt = sprintf('%s     (manual RR %.0f)', txt, hw.manual_rr); end
title(tl, txt, 'FontSize', 18, 'FontWeight', 'bold');
print(f, 'hardware_result.png', '-dpng', '-r200'); fprintf('Figure saved: hardware_result.png\n');

%% ======================= LOCAL FUNCTIONS =======================
function y = grid_(tx, vx, tg)                                       % same as extract_rmcs: 4 Hz grid, mean removed, z-scored
y = interp1(tx(:), vx(:), tg, 'pchip'); y = fillmissing(y, 'nearest'); y = y - mean(y); s = std(y); if s > 0, y = y/s; end
end
function y = zscore_(x)
x = x(:)' - mean(x); y = x/std(x);
end
function f = td_features_(x, epsfrac)                                % MAV WL ZC SSC SD AVP (same as rmc_classifier_selection.m)
x = x(:); dx = diff(x); e = epsfrac*std(x);
f = [mean(abs(x)), sum(abs(dx)), sum((x(1:end-1).*x(2:end) < 0) & (abs(dx) >= e)), sum((-dx(1:end-1).*dx(2:end)) > e^2), std(x), mean(x.^2)];
end
function [yp, sc] = gnb_(Xtr, ytr, Xte)                              % Gaussian naive Bayes, uniform priors; sc = P(component correct)
vfloor = 1e-6*var(Xtr, 0, 1) + 1e-12; ll = zeros(size(Xte, 1), 2);
for k = 0:1
    Xk = Xtr(ytr == k, :); mu = mean(Xk, 1); v = var(Xk, 0, 1) + vfloor;
    ll(:, k+1) = log(0.5) + sum(-0.5*log(2*pi*v) - (Xte - mu).^2 ./ (2*v), 2);
end
sc = 1 ./ (1 + exp(ll(:, 1) - ll(:, 2))); yp = double(ll(:, 2) > ll(:, 1));
end
