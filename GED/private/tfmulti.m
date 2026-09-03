function [pow, tout, freqs] = tfmulti(x, srate, opts)
% TFMULTI  Time-frequency power spectral density, by the multitaper method.
%
%   [psd, t, freqs] = tfmulti(x, srate)
%   [psd, t, freqs] = tfmulti(x, srate, Name, Value, ...)
%
%   Thomson's method (1982): taper the same stretch of signal several times over
%   with mutually orthogonal Slepian sequences, take a periodogram from each, and
%   average them. The tapers are the sequences most concentrated inside a
%   half-bandwidth W, so each periodogram leaks very little outside that band and
%   - being orthogonal - each is a nearly independent look at the same data.
%   Averaging K of them cuts the variance by about K without the smearing that
%   averaging neighbouring frequencies would cause.
%
%   How this differs from tfmorlet, which matters when reading the two side by
%   side:
%
%   Resolution is constant, not proportional. Every row is estimated with the
%   same bandwidth 2W = 2*tw/window Hz and the same window length, so 1 Hz and
%   30 Hz are measured alike. A wavelet bank instead holds the number of cycles
%   roughly fixed, so its resolution scales with frequency - sharp in time and
%   blunt in frequency at the top, the reverse at the bottom. Multitaper is the
%   better choice when the question is "which frequency", wavelets when it is
%   "exactly when".
%
%   The estimate is a power spectral density, in squared units of x per Hz,
%   where tfmorlet returns band power in squared units of x. The two therefore
%   differ by a scale factor that is not constant across frequency, and their
%   raw values are not interchangeable. Under the per-frequency normalisation
%   drawtfmap applies by default this makes no difference at all, since any
%   fixed factor divides out.
%
%   Variance is already low. K tapers make each estimate chi-squared with 2K
%   degrees of freedom rather than the 2 of a single wavelet estimate, which is
%   why this needs no equivalent of tfmorlet's smoothsigma.
%
%   Sub-windows, so nothing is skipped. Windows are placed at a hop of at most
%   half a window and averaged into the requested output bins. Sliding by a whole
%   output bin instead would, at the zoomed-out end where a bin spans minutes,
%   sample a few seconds and ignore everything between - the events in the gaps
%   would simply not appear.
%
%   Required
%   --------
%   x        Signal, a row vector. Single precision is fine.
%   srate    Sampling rate in Hz.
%
%   Name-Value
%   ----------
%   freqlim   [lo hi] of the frequency bank, in Hz. Default: [1 40]. Clipped to
%             the Nyquist frequency.
%   nfreq     Number of evenly spaced frequencies. Default: 40.
%   freqs     Explicit frequency vector, overriding freqlim and nfreq.
%   window    Length of the taper window, in seconds. Default: 4. This sets the
%             time resolution outright, and with tw the frequency resolution.
%   tw        Time-bandwidth product NW. Default: 3. The half-bandwidth is
%             tw/window Hz, so the resolution is 2*tw/window - 1.5 Hz at the
%             defaults. Raising it buys more usable tapers, and therefore less
%             variance, at the cost of a wider band.
%   ntapers   Number of tapers. 0 (default) takes 2*tw-1, the standard choice:
%             beyond that the sequences are no longer well concentrated and
%             start dragging in power from outside the band.
%   outrate   Output rate in Hz. Bins are floor(srate/outrate) samples wide, and
%             match tfmorlet's bins exactly, so the two put a column at the same
%             instant. 0 (default) keeps the input rate.
%   range     [i1 i2] sample indices to return, or a k-by-2 list of them,
%             concatenated. Context is taken from the real signal where there is
%             any, mirrored only at the ends of the recording.
%
%   Returns
%   -------
%   pow      nfreq-by-ntime single, one-sided PSD in units^2/Hz.
%   tout     Centre time of each column, in seconds from x(1).
%   freqs    Frequency of each row, in Hz.
%
%   Example
%   -------
%     [psd, t, f] = tfmulti(GED.comp(1, :), 250, 'window', 5, 'tw', 4);
%     imagesc(t, f, 10 * log10(psd ./ median(psd, 2))); axis xy
%
%   See also TFMORLET, DRAWTFMAP, PLOTGED.

arguments
    x            (1,:) {mustBeNumeric, mustBeReal}
    srate        (1,1) double {mustBePositive}
    opts.freqlim (1,2) double {mustBePositive} = [1 40]
    opts.nfreq   (1,1) double {mustBePositive} = 40
    opts.freqs         double = []
    opts.window  (1,1) double {mustBePositive} = 4
    opts.tw      (1,1) double {mustBePositive} = 3
    opts.ntapers (1,1) double {mustBeNonnegative} = 0
    opts.outrate (1,1) double {mustBeNonnegative} = 0
    opts.range         double = []
end

n = numel(x);

%% ---------------------------------------------------------------- the bank
freqs = opts.freqs;
if isempty(freqs)
    hi = min(opts.freqlim(2), srate / 2);
    lo = min(opts.freqlim(1), hi);
    freqs = linspace(lo, hi, opts.nfreq);
end
freqs = freqs(:).';
nfreq = numel(freqs);

%% ------------------------------------------------------------- the ranges
ranges = opts.range;
if isempty(ranges), ranges = [1 n]; end
ranges       = round(reshape(ranges, [], 2));
ranges(:, 1) = max(ranges(:, 1), 1);
ranges(:, 2) = min(ranges(:, 2), n);
ranges       = ranges(ranges(:, 2) >= ranges(:, 1), :);
if isempty(ranges)
    pow = zeros(nfreq, 0, 'single'); tout = zeros(1, 0);
    return
end
lens = ranges(:, 2) - ranges(:, 1) + 1;

%% ------------------------------------------------------- transform geometry
%%% Odd, so a window has a defined centre sample to hang the output column on.
nwin = max(3, round(opts.window * srate));
nwin = 2 * floor(nwin / 2) + 1;
offs = (-(nwin - 1) / 2 : (nwin - 1) / 2).';

K = opts.ntapers;
if K <= 0, K = max(1, round(2 * opts.tw - 1)); end
tapers = single(slepians(nwin, opts.tw, K));

dec = 1;
if opts.outrate > 0, dec = max(1, floor(srate / opts.outrate)); end

%%% Zero-pad the FFT until its bins are finer than the frequencies asked for,
%%% so the interpolation below reads a smooth curve rather than stepping between
%%% widely spaced bins. Padding adds no resolution - that is set by the window -
%%% only sampling of the same underlying estimate.
dfreq = srate / nwin;
if nfreq > 1, dfreq = min(dfreq, min(diff(freqs))); end
nfft  = 2^nextpow2(max(nwin, ceil(2 * srate / max(dfreq, eps))));
nbin  = nfft / 2 + 1;

%%% Linear interpolation from FFT bins onto the requested frequencies, as index
%%% and weight rather than a sparse matrix - sparse would force the whole
%%% spectrogram up to double on every multiply.
b  = freqs(:) * nfft / srate;
b0 = min(max(floor(b), 0), nbin - 2);
bw = single(b - b0);
i0 = b0 + 1;
i1 = b0 + 2;

%% ------------------------------------------------------------- the windows
nb   = ceil(lens / dec);
pow  = zeros(nfreq, sum(nb), 'single');
tout = zeros(1, sum(nb));

%%% Enough sub-windows that consecutive ones overlap by at least half, so every
%%% sample lands under the flat part of some taper. Capped, because past a point
%%% the extra looks are redundant and only cost time.
nsub = max(1, min(32, ceil(dec / max(1, floor(nwin / 2)))));

col = 0;
for r = 1:size(ranges, 1)
    o  = ranges(r, 1);
    nk = nb(r);

    %%% Sub-window centres, nsub of them spread evenly across each output bin.
    %%% Column-major order groups them by bin, which is what lets the averaging
    %%% below be a plain reshape.
    cen = round(o + (0:nk - 1) * dec + ((0:nsub - 1).' + 0.5) * (dec / nsub) - 0.5);
    cen = cen(:).';

    pw = zeros(nfreq, numel(cen), 'single');
    %%% Chunked, so the segment matrix and its transform stay bounded however
    %%% much of the recording was asked for.
    for c0 = 1:2048:numel(cen)
        cc  = c0:min(c0 + 2047, numel(cen));
        seg = single(x(reflectidx(cen(cc) + offs, n)));      % nwin x numel(cc)

        acc = zeros(nbin, numel(cc), 'single');
        for k = 1:K
            X   = fft(seg .* tapers(:, k), nfft, 1);
            X   = X(1:nbin, :);
            acc = acc + real(X).^2 + imag(X).^2;
        end
        %%% One-sided PSD: average over tapers, divide by the sampling rate, and
        %%% double to fold the negative frequencies in. The tapers carry unit
        %%% energy, so no window-power correction is needed on top.
        acc = acc * (2 / (K * srate));

        pw(:, cc) = (1 - bw) .* acc(i0, :) + bw .* acc(i1, :);
    end

    if nsub > 1
        pw = permute(mean(reshape(pw, nfreq, nsub, nk), 2), [1 3 2]);
    end

    pow(:, col + (1:nk)) = pw;
    tout(col + (1:nk))   = (o - 1 + (0:nk - 1) * dec + (dec - 1) / 2) / srate;
    col = col + nk;
end
end

% -------------------------------------------------------------------------
function E = slepians(N, NW, K)
% The first K discrete prolate spheroidal sequences, length N, half-bandwidth
% NW/N. Uses the Signal Processing Toolbox where it is there, and otherwise
% falls back - the diagnostic figures in this repo do not depend on a toolbox.
%
% The fallback is the standard construction (Percival & Walden 1993, ch. 8): the
% sequences are the eigenvectors of a symmetric tridiagonal matrix that commutes
% with the concentration kernel, which is far better conditioned than the kernel
% itself. Ordering by descending eigenvalue of that matrix gives them in order of
% descending concentration, which is the order wanted here.

if exist('dpss', 'file') == 2
    E = dpss(N, NW, K);
    return
end

i = (0:N - 1).';
d = ((N - 1 - 2 * i) / 2).^2 * cos(2 * pi * NW / N);
e = i(2:end) .* (N - i(2:end)) / 2;

[V, L] = eig(diag(d) + diag(e, 1) + diag(e, -1));
[~, ord] = sort(diag(L), 'descend');
E = V(:, ord(1:K));

%%% eig returns unit-norm columns already, so only the sign is left to fix. Even
%%% orders are made positive-going, odd orders positive at the start, which is
%%% the convention dpss uses - it costs nothing and keeps the two paths swappable.
for k = 1:K
    if mod(k - 1, 2) == 0
        if sum(E(:, k)) < 0, E(:, k) = -E(:, k); end
    elseif E(find(abs(E(:, k)) > eps, 1), k) < 0
        E(:, k) = -E(:, k);
    end
end
end
