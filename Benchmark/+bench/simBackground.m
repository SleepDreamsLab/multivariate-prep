function B = simBackground(fwd, stageEpoch, srate, epochLength, mode, opts)
% SIMBACKGROUND  Artifact-free background EEG with a known ground truth.
%
%   B = bench.simBackground(fwd, stageEpoch, srate, epochLength, 'synthetic', 'targetrms', rms)
%   B = bench.simBackground(fwd, stageEpoch, srate, epochLength, 'flat')
%
%   'synthetic'  nsources random cortical sources, each an independent 1/f^exponent noise
%                process (high-passed at 0.3 Hz), projected through the cortical leadfield,
%                plus white sensor noise. Every component lies in the span GEDAI keeps, so
%                anything GEDAI removes from it is over-cleaning by definition - the case the
%                raw (unpaired) SER/ARR variants are meant for.
%   'flat'       white sensor noise only. GEDAI's covariance-based thresholds degenerate on a
%                truly flat signal, so a small noise floor stands in for "flat".
%
%   Name-value
%   ----------
%   stages       stage digits that targetrms and exponent refer to        ([-2 -3 0 1])
%   targetrms    uV, per stage: median channel RMS after average reference; typically
%                measured on the real epochs so synthetic and real data share a scale.
%                [] leaves the unscaled leadfield amplitude                   ([])
%   exponent     1/f exponent per stage                                  ([2.0 2.4 1.5 1.2])
%   nsources     number of active cortical sources                        (2000)
%   sensornoise  uV RMS of white noise per channel                        (1 synthetic, 0.5 flat)
%
%   B is [nChan x nSamples] single, referenced to the recording reference.

arguments
    fwd struct
    stageEpoch (1,:) double
    srate (1,1) double {mustBePositive}
    epochLength (1,1) double {mustBePositive}
    mode {mustBeMember(mode, {'synthetic', 'flat'})}
    opts.stages (1,:) double = [-2 -3 0 1]
    opts.targetrms (1,:) double = []
    opts.exponent (1,:) double = [2.0 2.4 1.5 1.2]
    opts.nsources (1,1) double {mustBeInteger, mustBePositive} = 2000
    opts.sensornoise double = []
end

nChan        = size(fwd.Gctx, 1);
epochSamples = round(epochLength * srate);
nEpochs      = numel(stageEpoch);
B = zeros(nChan, nEpochs * epochSamples, 'single');

if strcmp(mode, 'flat')
    if isempty(opts.sensornoise), opts.sensornoise = 0.5; end
    B(:) = opts.sensornoise * randn(size(B), 'single');
    return
end
if isempty(opts.sensornoise), opts.sensornoise = 1; end

sources = randperm(size(fwd.Gctx, 2), min(opts.nsources, size(fwd.Gctx, 2)));
G = fwd.Gctx(:, sources);
G = G / sqrt(mean(G.^2, 'all'));

freqs = (0:epochSamples - 1) * srate / epochSamples;
freqs = min(freqs, srate - freqs);               % two-sided frequency of each FFT bin

for iStage = 1:numel(opts.stages)
    epochs = find(stageEpoch == opts.stages(iStage));
    if isempty(epochs), continue; end
    shaping = 1 ./ max(freqs, 0.3).^(opts.exponent(min(iStage, numel(opts.exponent))) / 2);
    shaping(freqs < 0.3) = 0;

    for iEpoch = epochs
        cols = (iEpoch - 1) * epochSamples + 1:iEpoch * epochSamples;
        sourceTs = real(ifft(fft(randn(numel(sources), epochSamples), [], 2) .* shaping, [], 2));
        B(:, cols) = single(G * sourceTs);
    end

    if ~isempty(opts.targetrms)
        cols = reshape((epochs - 1) * epochSamples + (1:epochSamples)', 1, []);
        current = median(sqrt(mean(double(bench.avgRef(B(:, cols))).^2, 2)));
        B(:, cols) = B(:, cols) * (opts.targetrms(iStage) / current);
    end
end
B = B + opts.sensornoise * randn(size(B), 'single');
end
