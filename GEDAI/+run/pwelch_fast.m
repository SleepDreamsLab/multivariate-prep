function [EpochPower, Frequencies, Times_h, KeepTime] = pwelch_fast(EEG, EpochLength, WindowLength, Overlap, KeepTime, MaxBytes)
% RUN.PWELCH_FAST  Batched Welch PSD per epoch — drop-in, faster replacement for run.run_pwelch.
%
%   [EpochPower, Frequencies, Times_h, KeepTime] = run.pwelch_fast(EEG)
%   [...] = run.pwelch_fast(EEG, EpochLength, WindowLength, Overlap)
%
%   Splits EEG.data into consecutive epochs and computes a Welch PSD per
%   epoch. Numerically identical to run.run_pwelch (which loops pwelch over
%   epochs via oscip.compute_power_on_epochs), but ~2x faster: all Welch
%   segments of an epoch block are gathered and transformed in a single
%   batched FFT instead of one pwelch call per epoch.
%
%   Inputs
%   ------
%   EEG          : EEGLAB struct; uses EEG.data (chans x samples) and EEG.srate.
%   EpochLength  : Epoch duration in seconds (default 30).
%   WindowLength : Welch window length in seconds (default 4). Also the nfft.
%   Overlap      : Welch overlap fraction, 0-1 (default 0.5).
%   KeepTime     : Optional timing struct to append to.
%   MaxBytes     : Working-buffer target per epoch block (default 256e6).
%                  Lower it if memory is tight; it only affects speed.
%
%   Outputs
%   -------
%   EpochPower   : chans x freqs x epochs, single. NOTE the dimension order —
%                  run.run_pwelch returned chans x epochs x freqs.
%   Frequencies  : 1 x freqs, single.
%   Times_h      : 1 x epochs, epoch start times in hours.
%   KeepTime     : KeepTime with .Pwelch set to the elapsed seconds.
%
%   NaN handling matches run.run_pwelch: a channel is NaN for a given epoch
%   iff it contains any NaN within that epoch.

arguments
    EEG
    EpochLength        = 30
    WindowLength       = 4
    Overlap            = 0.5
    KeepTime           = []
    MaxBytes     (1,1) = 256e6
end

fprintf('Computing power ...\n')
D = tic;

srate = EEG.srate;
Data  = EEG.data;

[nCh, nPts] = size(Data);
spe   = round(EpochLength * srate);           % samples per epoch
nEp   = floor(nPts / spe);
WP    = round(WindowLength * srate);          % window length, also nfft
step  = WP - round(WP * Overlap);
nSeg  = floor((spe - WP)/step) + 1;           % Welch segments fully inside an epoch
nF    = floor(WP/2) + 1;

% (0:nF-1)*srate/WP, not linspace(0,srate/2,nF): for odd nfft the last bin
% sits below Nyquist, and pwelch reports it that way.
Frequencies = single((0:nF-1) * (srate/WP));

win   = hanning(WP);
scale = 1/(srate * sum(win.^2));               % pwelch PSD normalisation
win   = single(win);

% Sample offsets of every Welch segment within one epoch (WP x nSeg)
off = (0:WP-1)' + (0:nSeg-1)*step + 1;

% Only pay for NaN bookkeeping if the recording actually has NaNs
hasNan = anynan(Data);

% chans x freqs x epochs — epoch-contiguous, which is also the output order
EpochPower = zeros(nCh, nF, nEp, 'single');

% Process epochs in blocks so the FFT working buffer stays near MaxBytes
bytesPerEp = WP * nSeg * nCh * 8;              % complex single
epChunk    = max(1, min(nEp, floor(MaxBytes / bytesPerEp)));

for ep0 = 1:epChunk:nEp
    ep1  = min(ep0 + epChunk - 1, nEp);
    nEpC = ep1 - ep0 + 1;

    X = Data(:, (ep0-1)*spe + 1 : ep1*spe).';   % (nEpC*spe) x nCh
    if ~isa(X, 'single'), X = single(X); end
    if hasNan
        nanMask = isnan(X);
        X(nanMask) = 0;                          % zero-fill so the FFT stays finite
    end

    % Gather every Welch segment: WP x nSeg x nEpC x nCh
    ridx = off + reshape((0:nEpC-1)*spe, 1, 1, nEpC);
    S    = reshape(X(ridx(:), :), WP, nSeg, nEpC, nCh);

    Y = fft(S .* win, WP, 1);

    % One-sided power, averaged over segments. Squaring the sliced real/imag
    % parts directly avoids materialising abs(Y) over the full spectrum.
    P = real(Y(1:nF,:,:,:)).^2 + imag(Y(1:nF,:,:,:)).^2;
    P = mean(P, 2) * scale;

    % Double everything but DC and, for even nfft, Nyquist
    if mod(WP, 2)
        P(2:end, :)   = P(2:end, :)   * 2;
    else
        P(2:end-1, :) = P(2:end-1, :) * 2;
    end

    EpochPower(:, :, ep0:ep1) = permute(reshape(P, nF, nEpC, nCh), [3 1 2]);

    if hasNan
        % A channel is NaN for an epoch iff it had any NaN anywhere in it
        bad = reshape(any(reshape(nanMask.', nCh, spe, nEpC), 2), nCh, nEpC);
        blk = EpochPower(:, :, ep0:ep1);
        blk(repmat(reshape(bad, nCh, 1, nEpC), 1, nF, 1)) = NaN;
        EpochPower(:, :, ep0:ep1) = blk;
    end
end

KeepTime.Pwelch = toc(D);

%%% Times vector
Times_h = (1:nEp) .* EpochLength/60/60;

fprintf('Computing power took %.0fs!\n', KeepTime.Pwelch)

end
