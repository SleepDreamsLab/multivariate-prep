function [det, thresholds] = detectSlowWaves(x, srate, mask, opts)
% DETECTSLOWWAVES  Zero-crossing slow-wave / K-complex detector with optionally frozen thresholds.
%
%   [det, thresholds] = bench.detectSlowWaves(x, srate, mask)
%   det = bench.detectSlowWaves(x, srate, mask, 'thresholds', thresholds)
%
%   Zero-crossing detection after Massimini et al. (2004, J Neurosci 24:6862): in the
%   band-passed signal, a candidate is a negative half-wave (down zero crossing to up zero
%   crossing) of plausible length followed by its positive half-wave up to the next down
%   crossing. Amplitude criteria are relative to the candidate distribution, as in many
%   later adaptations (e.g. Staresina et al. 2015, Nat Neurosci 18:1679, who kept the
%   largest quarter): a candidate is kept when its negative peak is at least negfactor x
%   the mean candidate negative peak and its peak-to-peak at least p2pfactor x the mean
%   candidate peak-to-peak. Absolute AASM-style values would not transfer to
%   average-referenced 256-channel data.
%
%   As for bench.detectSpindles, bidsfun_gedai_benchmark calibrates the thresholds once on
%   the uncleaned input and freezes them for the cleaned signals, so amplitude loss shows up
%   as missed events.
%
%   Inputs
%   ------
%   x          [1 x n] uV
%   mask       [1 x n] logical; a candidate must lie entirely inside it        (true(1,n))
%
%   Name-value
%   ----------
%   band          Hz                                                      ([0.5 4])
%   negduration   s, [min max] length of the negative half-wave            ([0.25 1.0])
%   maxduration   s, maximum length of the whole wave                      (2.5)
%   negfactor     multiple of the mean candidate negative peak             (1)
%   p2pfactor     multiple of the mean candidate peak-to-peak               (1)
%   thresholds    [negPeak p2p] uV (negPeak negative); [] calibrates on x  ([])
%
%   Output
%   ------
%   det         struct of column vectors: onset, offset (samples), duration (s), negpeak
%               (uV), negpeaksample, amplitude (peak-to-peak uV)
%   thresholds  [negPeak p2p] used

arguments
    x (1,:) double
    srate (1,1) double {mustBePositive}
    mask (1,:) logical = true(size(x))
    opts.band (1,2) double = [0.5 4]
    opts.negduration (1,2) double = [0.25 1.0]
    opts.maxduration (1,1) double = 2.5
    opts.negfactor (1,1) double = 1
    opts.p2pfactor (1,1) double = 1
    opts.thresholds double = []
end

filtered = bandpassZeroPhase(x, srate, opts.band, 2);
positive = filtered >= 0;
downZc = find(positive(1:end-1) & ~positive(2:end));     % last sample before going negative
upZc   = find(~positive(1:end-1) & positive(2:end));

%%% For each down crossing, the next up crossing and the down crossing after that
allDown = downZc;
nextUp = discretize(downZc, [-Inf, upZc, Inf]);          % bin b: upZc(b) is the next up crossing
valid  = nextUp <= numel(upZc);
downZc = downZc(valid);
up     = upZc(nextUp(valid));
nextDown = discretize(up, [-Inf, allDown, Inf]);         % bin b: allDown(b) is the next down crossing
valid  = nextDown <= numel(allDown);
downZc = downZc(valid);
up     = up(valid);
down2  = allDown(nextDown(valid));

negLength   = (up - downZc) / srate;
totalLength = (down2 - downZc) / srate;
outside     = cumsum(~mask);
inMask      = mask(downZc) & (outside(down2) - outside(downZc)) == 0;
isCandidate = negLength >= opts.negduration(1) & negLength <= opts.negduration(2) & ...
    totalLength <= opts.maxduration & inMask;
downZc = downZc(isCandidate);
up     = up(isCandidate);
down2  = down2(isCandidate);

n = numel(downZc);
negPeak = zeros(n, 1); negSample = zeros(n, 1); p2p = zeros(n, 1);
for k = 1:n
    [negPeak(k), iMin] = min(filtered(downZc(k):up(k)));
    negSample(k) = downZc(k) + iMin - 1;
    p2p(k) = max(filtered(up(k):down2(k))) - negPeak(k);
end

thresholds = opts.thresholds;
if isempty(thresholds)
    if n == 0
        thresholds = [-Inf Inf];
    else
        thresholds = [opts.negfactor * mean(negPeak), opts.p2pfactor * mean(p2p)];
    end
end
keep = negPeak <= thresholds(1) & p2p >= thresholds(2);

det = struct('onset', downZc(keep)' + 1, 'offset', down2(keep)', ...
    'duration', (down2(keep)' - downZc(keep)') / srate, 'negpeak', negPeak(keep), ...
    'negpeaksample', negSample(keep), 'amplitude', p2p(keep));
end
