function [pow, tout, freqs] = tfmorlet(x, srate, opts)
% TFMORLET  Time-frequency power, by complex Morlet wavelet convolution.
%
%   [pow, t, freqs] = tfmorlet(x, srate)
%   [pow, t, freqs] = tfmorlet(x, srate, Name, Value, ...)
%
%   The textbook transform (Cohen 2014, ch. 12-13): convolve the signal with a
%   bank of complex Morlet wavelets and take the squared magnitude. What is worth
%   knowing is how the practical decisions are made here.
%
%   Evenly spaced frequencies, so the bank can be drawn on a linear frequency
%   axis - image display assumes evenly spaced rows, and the spacing of the bank
%   and the scale of the axis are therefore not independent choices. Note what
%   this costs: a wavelet's bandwidth grows with its centre frequency, so an
%   evenly spaced bank oversamples the top of the range, where neighbouring rows
%   overlap heavily and carry nearly the same estimate. That is redundancy, not
%   error. Pass freqs explicitly for a log-spaced bank, which matches the
%   spacing to the resolution but needs a log axis to be drawn honestly.
%
%   A cycle count that rises with frequency, linearly across the bank. Fixing it
%   low buys time precision the low frequencies cannot use and leaves the high
%   ones smeared across bands; fixing it high makes the low frequencies span
%   tens of seconds. Rising from 3 to 10 is the usual compromise.
%
%   Unit peak gain. Each wavelet is divided by the maximum of its own spectrum,
%   so power comes out in the squared units of x rather than scaling with the
%   wavelet's length - which is what makes rows at different frequencies (and
%   different cycle counts) comparable at all.
%
%   Reflected edges. The signal is mirrored where a wavelet runs off its end,
%   so it never convolves against a step down to zero.
%
%   FFT convolution in overlap-save blocks. Doing it in one transform over a
%   whole night would need an FFT of millions of points per frequency; blocks
%   bound the memory to one block by nfreq and let the wavelet spectra be built
%   once and reused for every one of them.
%
%   Required
%   --------
%   x        Signal, a row vector. Single precision is fine - it is converted
%            block by block, never in bulk.
%   srate    Sampling rate in Hz.
%
%   Name-Value
%   ----------
%   freqlim   [lo hi] of the frequency bank, in Hz. Default: [1 40]. Clipped to
%             the Nyquist frequency.
%   nfreq     Number of evenly spaced frequencies. Default: 40.
%   freqs     Explicit frequency vector, overriding freqlim and nfreq.
%   cycles    [lo hi] wavelet cycles, spread linearly across the bank.
%             Default: [3 10].
%   outrate   Output rate in Hz. The power is averaged into bins of
%             floor(srate/outrate) samples - averaged, not sampled, so the bins
%             are their own anti-alias filter. 0 (default) keeps the input rate.
%   smoothsigma  Smooth each row in time by this many of its own wavelet widths
%             (sigma). A single power estimate is chi-squared with 2 degrees of
%             freedom - it scatters over some 25 dB whatever the signal does -
%             and displaying that raw spends the whole colour scale on the
%             scatter. Smoothing to one sigma throws away no detail the wavelet
%             could have resolved in the first place, since the row's power
%             cannot change faster than the wavelet that measured it. 0
%             (default) returns the unsmoothed power. Whatever binning outrate
%             already did counts towards the total, so the result carries the
%             same smoothing at every output rate.
%   range     [i1 i2] sample indices to return, or a k-by-2 list of them, which
%             comes back as one concatenated result. Context on either side is
%             taken from the real signal wherever there is any, so a window
%             computed on its own matches what the full-length transform would
%             have put there. Default: the whole signal.
%
%   Returns
%   -------
%   pow      nfreq-by-ntime single, |analytic signal|^2.
%   tout     Centre time of each column, in seconds from x(1). Not contiguous
%            when several ranges were asked for.
%   freqs    Frequency of each row, in Hz.
%
%   Example
%   -------
%     [pow, t, f] = tfmorlet(GED.comp(1, :), 250, 'freqlim', [1 25], 'outrate', 20);
%     imagesc(t, 1:numel(f), 10 * log10(pow ./ median(pow, 2)));
%
%   See also DRAWTFMAP, PLOTGED.

arguments
    x             (1,:) {mustBeNumeric, mustBeReal}
    srate         (1,1) double {mustBePositive}
    opts.freqlim  (1,2) double {mustBePositive} = [1 40]
    opts.nfreq    (1,1) double {mustBePositive} = 40
    opts.freqs          double = []
    opts.cycles   (1,2) double {mustBePositive} = [3 10]
    opts.outrate  (1,1) double {mustBeNonnegative} = 0
    opts.range          double = []
    opts.smoothsigma (1,1) double {mustBeNonnegative} = 0
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

if nfreq > 1
    ncyc = linspace(opts.cycles(1), opts.cycles(2), nfreq);
else
    ncyc = mean(opts.cycles);
end
sigma = ncyc ./ (2 * pi * freqs);            % Gaussian width, in seconds

%% ------------------------------------------------------------- the ranges
ranges = opts.range;
if isempty(ranges), ranges = [1 n]; end
ranges      = round(reshape(ranges, [], 2));
ranges(:, 1) = max(ranges(:, 1), 1);
ranges(:, 2) = min(ranges(:, 2), n);
ranges      = ranges(ranges(:, 2) >= ranges(:, 1), :);
if isempty(ranges)
    pow = zeros(nfreq, 0, 'single'); tout = zeros(1, 0);
    return
end
lens = ranges(:, 2) - ranges(:, 1) + 1;

%% ------------------------------------------------------- transform geometry
%%% 3.5 sigma of the widest wavelet: past that the Gaussian is down by a factor
%%% of ~2000 and its tail sits below the numerical noise of the convolution.
halfwave = ceil(3.5 * max(sigma) * srate);
nWave    = 2 * halfwave + 1;

dec = 1;
if opts.outrate > 0, dec = max(1, floor(srate / opts.outrate)); end

%%% Each block costs nWave-1 samples of overlap whatever its size, so short
%%% blocks spend most of their FFT on padding; large ones buy nothing once the
%%% overlap is amortised, and cost memory. A few times the wavelet, and never
%%% more output than was actually asked for.
nWant = min(max(lens), max(4 * nWave, 2^14));
nConv = 2^nextpow2(nWant + nWave - 1);
nData = floor((nConv - nWave + 1) / dec) * dec;    % whole output bins only
if nData < dec
    nConv = 2^nextpow2(nWave + dec);
    nData = floor((nConv - nWave + 1) / dec) * dec;
end

%%% The wavelet spectra, built once here and reused by every block below - the
%%% whole reason for keeping one FFT length for the entire bank.
wavtime = (-halfwave:halfwave) / srate;
waveX   = complex(zeros(nConv, nfreq));
for fi = 1:nfreq
    cmw = exp(2i * pi * freqs(fi) * wavtime - wavtime.^2 / (2 * sigma(fi)^2));
    %%% Divided by the analytic peak of the Gaussian's spectrum, srate*sigma*
    %%% sqrt(2*pi), rather than by max(abs(fft(...))). The sampled maximum lands
    %%% wherever the FFT grid happens to fall next to the peak, which makes the
    %%% scaling depend on the block length - so the same stretch of signal would
    %%% come out a fraction of a percent different depending on how much of the
    %%% recording was asked for. The analytic value does not move.
    waveX(:, fi) = fft(cmw, nConv).' / (srate * sigma(fi) * sqrt(2 * pi));
end
%%% Built in double and demoted here: the wavelets themselves want the precision,
%%% the convolution does not. A single-precision FFT is roughly twice as fast and
%%% moves half the memory, and the result - a power estimate that scatters over
%%% decibels - is nowhere near single's 1e-7.
waveX = complex(single(real(waveX)), single(imag(waveX)));

%% -------------------------------------------------------------- convolution
nb   = ceil(lens / dec);
pow  = zeros(nfreq, sum(nb), 'single');
tout = zeros(1, sum(nb));

%%% Overlap-save: a block reads nConv samples and yields nData outputs that are
%%% free of the FFT's circular wrap. Only samples from nWave on are untouched by
%%% the wrap, and the wavelet is centred in its own array, so keep(1) is the
%%% filtered value at the block's first output sample.
keep = (nWave - 1) + (1:nData);
col  = 0;
for r = 1:size(ranges, 1)
    o    = ranges(r, 1);
    left = lens(r);
    while left > 0
        L   = min(nData, left);
        seg = single(x(reflectidx(o - halfwave + (0:nConv - 1), n)));
        A   = ifft(fft(seg(:)) .* waveX, [], 1);
        A   = A(keep(1:L), :);
        %%% real^2 + imag^2 rather than abs(.)^2, which would take a square root
        %%% only to square it again - and this is the inner loop.
        P   = real(A).^2 + imag(A).^2;                  % L x nfreq

        k = ceil(L / dec);
        if dec == 1
            Pb = P.';
        else
            %%% NaN pads the last, partial bin rather than zeros, which would
            %%% pull its average down towards nothing.
            P(end + 1:k * dec, :) = NaN;
            Pb = permute(mean(reshape(P, dec, k, nfreq), 1, 'omitnan'), [3 2 1]);
        end
        pow(:, col + (1:k)) = Pb;
        tout(col + (1:k))   = (o - 1 + (0:k - 1) * dec + (dec - 1) / 2) / srate;

        col  = col + k;
        o    = o + L;
        left = left - L;
    end
end

%% ---------------------------------------------------------------- smoothing
%%% After the binning, not before it: the bins have already done part of the
%%% job, so only the shortfall is left to do here, and the total comes out the
%%% same whatever outrate was asked for. Doing it here also keeps it away from
%%% the block seams - the bins are contiguous across them, the blocks are not.
if opts.smoothsigma > 0
    w = max(1, round(opts.smoothsigma * sigma * srate / dec));
    if any(w > 1)
        stop  = cumsum(nb(:).');
        start = [1 stop(1:end - 1) + 1];
        for r = 1:numel(start)
            %%% One range at a time: they are separate stretches of recording
            %%% and smoothing across the join would mix them.
            cols = start(r):stop(r);
            for fi = find(w > 1)
                pow(fi, cols) = movmean(pow(fi, cols), w(fi));
            end
        end
    end
end
end
