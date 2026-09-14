function sel = selectEpochs(fdtFile, nbchan, pnts, srate, epochLength, scoring, opts)
% SELECTEPOCHS  Choose the benchmark epochs per sleep stage.
%
%   sel = bench.selectEpochs(fdtFile, nbchan, pnts, srate, epochLength, scoring)
%   sel = bench.selectEpochs(..., 'stages', [-2 -3 0 1], 'nepochs', 40, 'selectmode', 'clean')
%
%   The benchmark injects into a subset of each stage rather than the whole night: 40
%   epochs x 4 stages x 256 channels is ~1.2 GB per copy in single, and the benchmark holds
%   about six copies (background, baseline GEDAI output, truth, contaminated input, GEDAI
%   output, scoring buffers). Peak memory therefore scales with nepochs; see
%   bidsfun_gedai_benchmark.
%
%   Inputs
%   ------
%   fdtFile      EEGLAB .fdt of the input recording (float32, multiplexed)
%   nbchan, pnts header of that recording
%   scoring      [1 x nEpochs] stage digits (scoreloader convention)
%
%   Name-value
%   ----------
%   stages              stage digits to draw from                          ([-2 -3 0 1])
%   nepochs             epochs per stage; fewer if the stage is shorter      (40)
%   selectmode          'clean'  - read a random candidate pool and keep the epochs with
%                                  the lowest artifact score (below)
%                       'random' - random epochs, nothing read
%                       'first'  - the first epochs of the stage             ('clean')
%   candidatepool       'clean' reads candidatepool * nepochs candidates per stage, which
%                       bounds the I/O on a network share                     (4)
%   excludetransitions  skip epochs whose neighbours are scored differently - stage
%                       boundaries are where arousals and movement live      (true)
%
%   Artifact score (selectmode 'clean'): per epoch and channel, log variance, log variance
%   of the first difference (a high-frequency/EMG proxy) and log max |x| after demeaning;
%   each robustly z-scored per channel across the pool (median/MAD); the score is the 90th
%   percentile across channels of the largest of the three z-values. Low is clean.
%
%   Output
%   ------
%   sel.epochs     [1 x n] selected epoch indices in the recording, chronological
%   sel.stage      [1 x n] their stage digits
%   sel.score      [1 x n] artifact score (NaN unless selectmode 'clean')
%   sel.available  [1 x numel(stages)] candidate epochs per stage after exclusions

arguments
    fdtFile char
    nbchan (1,1) double
    pnts (1,1) double
    srate (1,1) double
    epochLength (1,1) double
    scoring (1,:) double
    opts.stages (1,:) double = [-2 -3 0 1]
    opts.nepochs (1,1) double {mustBeInteger, mustBePositive} = 40
    opts.selectmode {mustBeMember(opts.selectmode, {'clean', 'random', 'first'})} = 'clean'
    opts.candidatepool (1,1) double {mustBePositive} = 4
    opts.excludetransitions (1,1) logical = true
end

epochSamples = round(epochLength * srate);
nEpochsFile  = min(numel(scoring), floor(pnts / epochSamples));
scoring      = scoring(1:nEpochsFile);

epochs = []; stages = []; scores = [];
sel.available = zeros(1, numel(opts.stages));

fid = -1;
if strcmp(opts.selectmode, 'clean')
    fid = fopen(fdtFile, 'r');
    if fid < 0
        error('bench:selectEpochs:open', 'Cannot open %s.', fdtFile);
    end
end
closer = onCleanup(@() closeIfOpen(fid));

for iStage = 1:numel(opts.stages)
    candidates = find(scoring == opts.stages(iStage));
    if opts.excludetransitions && ~isempty(candidates)
        previous = [NaN, scoring(1:end-1)];
        next     = [scoring(2:end), NaN];
        stable   = (previous(candidates) == opts.stages(iStage) | isnan(previous(candidates))) & ...
                   (next(candidates) == opts.stages(iStage) | isnan(next(candidates)));
        if any(stable)
            candidates = candidates(stable);
        end
    end
    sel.available(iStage) = numel(candidates);
    n = min(opts.nepochs, numel(candidates));
    if n == 0
        continue
    end

    switch opts.selectmode
        case 'first'
            pick = candidates(1:n);
            pickScore = nan(1, n);
        case 'random'
            pick = candidates(sort(randperm(numel(candidates), n)));
            pickScore = nan(1, n);
        case 'clean'
            pool = candidates(randperm(numel(candidates), min(numel(candidates), ceil(opts.candidatepool * n))));
            poolScore = artifactScore(fid, nbchan, pool, epochSamples);
            [~, order] = sort(poolScore);
            pick = pool(order(1:n));
            pickScore = poolScore(order(1:n));
        otherwise
            error('bench:selectEpochs:mode', 'Unknown selectmode %s.', opts.selectmode);
    end
    epochs = [epochs, pick(:)'];                             %#ok<AGROW>
    stages = [stages, repmat(opts.stages(iStage), 1, n)];    %#ok<AGROW>
    scores = [scores, pickScore(:)'];                        %#ok<AGROW>
end

[sel.epochs, order] = sort(epochs);
sel.stage = stages(order);
sel.score = scores(order);
end

% -------------------------------------------------------------------------
function score = artifactScore(fid, nbchan, pool, epochSamples)
features = zeros(nbchan, numel(pool), 3);
for k = 1:numel(pool)
    x = double(readOneEpoch(fid, nbchan, pool(k), epochSamples));
    x = x - mean(x, 2);
    features(:, k, 1) = log(var(x, 0, 2));
    features(:, k, 2) = log(var(diff(x, 1, 2), 0, 2));
    features(:, k, 3) = log(max(abs(x), [], 2));
end
features(~isfinite(features)) = NaN;                 % flat channels carry no information
centre = median(features, 2, 'omitnan');
spread = 1.4826 * median(abs(features - centre), 2, 'omitnan') + eps;
z = max((features - centre) ./ spread, [], 3);
score = prctile(z, 90, 1);
end

% -------------------------------------------------------------------------
function closeIfOpen(fid)
if fid >= 0
    fclose(fid);
end
end
