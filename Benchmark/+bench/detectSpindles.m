function [det, threshold] = detectSpindles(x, srate, mask, opts)
% DETECTSPINDLES  RMS-threshold spindle detector with an optionally frozen threshold.
%
%   [det, threshold] = bench.detectSpindles(x, srate, mask)
%   det = bench.detectSpindles(x, srate, mask, 'threshold', threshold)
%
%   The detector of Mölle et al. (2002, J Neurosci 22:10941): band-pass, moving RMS over
%   200 ms, threshold at mean + 1.5 SD of the RMS within mask, events 0.5-3 s long.
%
%   Why a frozen threshold: a threshold calibrated on each signal separately is relative to
%   that signal, so a cleaning step that shrinks everything by the same factor would leave
%   the detections unchanged and hide exactly the amplitude loss the benchmark measures.
%   bidsfun_gedai_benchmark calibrates once on the uncleaned input and applies that value
%   to the GEDAI output and to both baseline signals.
%
%   Inputs
%   ------
%   x          [1 x n] uV
%   mask       [1 x n] logical, samples eligible for calibration and detection (true(1,n))
%
%   Name-value
%   ----------
%   band         Hz                                                   ([11 16])
%   rmswindow    s                                                    (0.2)
%   thresholdsd  SDs above the mean RMS                               (1.5)
%   duration     s, [min max]                                         ([0.5 3])
%   threshold    uV RMS; [] calibrates on x                           ([])
%
%   Output
%   ------
%   det        struct of column vectors: onset, offset (samples), duration (s),
%              amplitude (peak-to-peak uV of the band-passed signal), freq (Hz, from zero
%              crossings)
%   threshold  the RMS threshold used (uV)

arguments
    x (1,:) double
    srate (1,1) double {mustBePositive}
    mask (1,:) logical = true(size(x))
    opts.band (1,2) double = [11 16]
    opts.rmswindow (1,1) double {mustBePositive} = 0.2
    opts.thresholdsd (1,1) double = 1.5
    opts.duration (1,2) double = [0.5 3]
    opts.threshold double = []
end

filtered = bandpassZeroPhase(x, srate, opts.band, 4);
rmsEnv = sqrt(movmean(filtered.^2, max(1, round(opts.rmswindow * srate))));

threshold = opts.threshold;
if isempty(threshold)
    values = rmsEnv(mask);
    threshold = mean(values) + opts.thresholdsd * std(values);
end

above  = rmsEnv > threshold & mask;
edges  = diff([false, above, false]);
starts = find(edges == 1);
ends   = find(edges == -1) - 1;
len    = (ends - starts + 1) / srate;
keep   = len >= opts.duration(1) & len <= opts.duration(2);
starts = starts(keep);
ends   = ends(keep);

n = numel(starts);
det = struct('onset', starts(:), 'offset', ends(:), 'duration', (ends(:) - starts(:) + 1) / srate, ...
    'amplitude', zeros(n, 1), 'freq', zeros(n, 1));
for k = 1:n
    segment = filtered(starts(k):ends(k));
    det.amplitude(k) = max(segment) - min(segment);
    det.freq(k) = nnz(diff(sign(segment)) ~= 0) / 2 / det.duration(k);
end
end
