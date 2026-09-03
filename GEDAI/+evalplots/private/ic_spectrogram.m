function [Pwr, Freqs, Tsec] = ic_spectrogram(x, srate, opts)
% IC_SPECTROGRAM  Welch spectrogram of one component's activation.
%
%   [Pwr, Freqs, Tsec] = ic_spectrogram(x, srate)
%   [Pwr, Freqs, Tsec] = ic_spectrogram(x, srate, Name, Value, ...)
%
%   The same estimator run.pwelch_fast uses - Hann-tapered, one-sided PSD in
%   units^2/Hz - but sliding over the recording rather than averaged within
%   fixed epochs, so the result is a time-frequency map instead of one spectrum
%   per epoch. Written here rather than reusing run.pwelch_fast because that one
%   collapses every window inside an epoch into a single column, which is
%   exactly the time resolution this plot needs to keep.
%
%   Two things are worth knowing about how it is done:
%
%   Windows are averaged into output columns after the transform, not sampled.
%   A whole night at 2 s hops is some 14000 windows and no figure tile can draw
%   that many; taking every n-th window instead would show a few seconds out of
%   every minute and silently drop everything between - which for artefact
%   components, the whole point of the plot, is the part that matters.
%
%   NaN columns stay NaN. A window containing any NaN is zeroed before the FFT
%   so the transform stays finite, and its power is then set to NaN, matching
%   run.pwelch_fast's per-epoch rule.
%
%   Required
%   --------
%   x        Component activation, a row vector.
%   srate    Sampling rate in Hz.
%
%   Name-Value
%   ----------
%   Window    Taper window length in seconds. Default: 4 - the same as the
%             evaluation figures' Welch window, so the frequency resolution
%             matches what the PSD panels show.
%   Overlap   Window overlap fraction, 0-1. Default: 0.5.
%   FreqLim   [lo hi] frequency range kept, in Hz. Default: [0.5 40].
%   MaxCols   Upper bound on output columns. Default: 600.
%
%   Returns
%   -------
%   Pwr      nFreq-by-nCol single, PSD in units^2/Hz.
%   Freqs    1-by-nFreq frequency vector, Hz.
%   Tsec     1-by-nCol centre time of each column, in seconds from x(1).

arguments
    x            (1,:) {mustBeNumeric, mustBeReal}
    srate        (1,1) double {mustBePositive}
    opts.Window  (1,1) double {mustBePositive} = 4
    opts.Overlap (1,1) double {mustBeNonnegative, mustBeLessThan(opts.Overlap, 1)} = 0.5
    opts.FreqLim (1,2) double = [0.5 40]
    opts.MaxCols (1,1) double {mustBePositive} = 600
end

x = single(x);
n = numel(x);

WP   = max(8, round(opts.Window * srate));
step = max(1, WP - round(WP * opts.Overlap));
nSeg = floor((n - WP) / step) + 1;
if nSeg < 1
    Pwr = zeros(0, 0, 'single'); Freqs = zeros(1, 0, 'single'); Tsec = zeros(1, 0);
    return
end

%%% (0:nF-1)*srate/WP, not linspace(0, srate/2, nF): for an odd window the last
%%% bin sits below Nyquist, and pwelch reports it that way.
nF   = floor(WP / 2) + 1;
fAll = single((0:nF-1) * (srate / WP));
keep = fAll >= opts.FreqLim(1) & fAll <= opts.FreqLim(2);
if ~any(keep), keep(:) = true; end
Freqs = fAll(keep);

win   = hanning(WP);
scale = 1 / (srate * sum(win.^2));
win   = single(win);

P = zeros(nnz(keep), nSeg, 'single');

%%% Segments in blocks, so the gathered window matrix and its transform stay
%%% bounded however long the recording is.
segBlock = max(1, floor(4e6 / WP));
for s0 = 1:segBlock:nSeg
    s1  = min(s0 + segBlock - 1, nSeg);
    idx = (0:WP-1)' + ((s0:s1) - 1) * step + 1;
    S   = x(idx);                                  % WP x nBlock

    bad = any(isnan(S), 1);
    if any(bad), S(:, bad) = 0; end

    Y  = fft(S .* win, WP, 1);
    Ps = (real(Y(1:nF, :)).^2 + imag(Y(1:nF, :)).^2) * scale;

    %%% Double everything but DC and, for an even window, Nyquist.
    if mod(WP, 2)
        Ps(2:end, :)   = Ps(2:end, :)   * 2;
    else
        Ps(2:end-1, :) = Ps(2:end-1, :) * 2;
    end

    Ps = Ps(keep, :);
    if any(bad), Ps(:, bad) = NaN; end
    P(:, s0:s1) = Ps;
end

tSeg = ((0:nSeg-1) * step + (WP - 1) / 2) / srate;

%%% Average neighbouring windows down to at most MaxCols columns. NaN pads the
%%% last, partial bin rather than zeros, which would pull its mean towards
%%% nothing rather than leaving it out.
binw = max(1, ceil(nSeg / opts.MaxCols));
if binw > 1
    nk  = ceil(nSeg / binw);
    pad = nk * binw - nSeg;
    P(:, end+1:end+pad)  = NaN;
    tSeg(end+1:end+pad)  = NaN;
    Pwr  = permute(mean(reshape(P, size(P, 1), binw, nk), 2, 'omitnan'), [1 3 2]);
    Tsec = mean(reshape(tSeg, binw, nk), 1, 'omitnan');
else
    Pwr  = P;
    Tsec = tSeg;
end
end
