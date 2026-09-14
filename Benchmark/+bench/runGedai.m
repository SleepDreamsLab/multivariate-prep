function Xhat = runGedai(EEG, data, stageEpoch, groups, r, lfCOV, touched, baseline, opts)
% RUNGEDAI  Run GEDAI stage-group-wise on benchmark data, re-running only groups that changed.
%
%   Bhat = bench.runGedai(EEG, B, stageEpoch, groups, r, lfCOV, true(1, nEpochs), single([]))
%   Xhat = bench.runGedai(EEG, X, stageEpoch, groups, r, lfCOV, touched, Bhat)
%
%   GEDAI runs through run.GEDAI_StageSpecific with the same arguments bidsfun_gedai passes,
%   so the benchmark exercises the production code path. A stage group is cleaned in one
%   pass with group-level thresholds, so a group is re-run whole whenever any of its epochs
%   was touched by an injection; groups that nothing touched are identical to the baseline
%   run and are copied from it instead of being cleaned again.
%
%   Inputs
%   ------
%   EEG         EEGLAB header of the benchmark recording (chanlocs, srate); data is ignored
%   data        [nChan x nEpochs*epochSamples] single, average-referenced
%   stageEpoch  [1 x nEpochs] stage digits
%   groups      .stageLogic (cell of stage digit vectors), .mode and .modeBB (cells of
%               GEDAI mode strings), as built by bidsfun_gedai_benchmark
%   r           GEDAI run configuration (gedai.defaultRuns entry)
%   lfCOV       leadfield covariance (gedai.loadrefcov)
%   touched     [1 x nEpochs] logical, epochs that differ from the baseline data
%   baseline    GEDAI output for the baseline data; single([]) runs every group
%
%   Name-value
%   ----------
%   epochlength  s                                                          (30)
%   pooltype     'Processes' or 'Threads', passed to GEDAI_StageSpecific    ('Processes')
%
%   Outputs
%   -------
%   Xhat        [nChan x nSamples] single

arguments
    EEG struct
    data single
    stageEpoch (1,:) double
    groups struct
    r struct
    lfCOV double
    touched (1,:) logical
    baseline single = single([])
    opts.epochlength (1,1) double = 30
    opts.pooltype {mustBeMember(opts.pooltype, {'Processes', 'Threads'})} = 'Processes'
end

epochSamples = round(opts.epochlength * EEG.srate);
present = cellfun(@(s) any(ismember(stageEpoch, s)), groups.stageLogic);
if isempty(baseline)
    runGroup = present;
    Xhat = zeros(size(data), 'single');
else
    runGroup = present & cellfun(@(s) any(touched & ismember(stageEpoch, s)), groups.stageLogic);
    Xhat = baseline;
end
if ~any(runGroup)
    return
end

epochs = find(ismember(stageEpoch, [groups.stageLogic{runGroup}]));
cols = reshape((epochs - 1) * epochSamples + (1:epochSamples)', 1, []);

EEGsub = EEG;
EEGsub.data     = data(:, cols);
EEGsub.nbchan   = size(EEGsub.data, 1);
EEGsub.pnts     = size(EEGsub.data, 2);
EEGsub.trials   = 1;
EEGsub.xmin     = 0;
EEGsub.xmax     = (EEGsub.pnts - 1) / EEGsub.srate;
EEGsub.times    = [];
EEGsub.event    = [];
EEGsub.urevent  = [];
EEGsub.epoch    = [];

nGroups = nnz(runGroup);
%%% Only the cleaned output: with nargout >= 2 GEDAI_StageSpecific also assembles the
%%% uncleaned data, a full extra copy the benchmark already holds.
EEGclean = run.GEDAI_StageSpecific(EEGsub, stageEpoch(epochs), ...
    groups.stageLogic(runGroup), struct(), ...
    'EpochLength',                opts.epochlength, ...
    'GEDAIMode',                  groups.mode(runGroup), ...
    'GEDAIModeBB',                groups.modeBB(runGroup), ...
    'GEDAIEpochSize',             r.GEDAIEpochSize, ...
    'GEDAILowCutOffFreq',         r.GEDAILowCutOffFreq, ...
    'BBEpochSize',                r.GEDAIBroadbandEpochSize, ...
    'BroadbandOnly',              r.broadbandOnly, ...
    'GEDAIEnovaChannelThreshold', r.GEDAIEnovaChannelThreshold, ...
    'PercentileThreshold',        r.percentileThreshold, ...
    'BBMinThreshold',             r.BBMinThreshold, ...
    'ComputeSENSAI',              r.computeSENSAI, ...
    'ICAtype',                    'none', ...
    'PoolType',                   opts.pooltype, ...
    'RefCOV',                     repmat({lfCOV}, 1, nGroups));

if ~isequal(size(EEGclean.data), [size(data, 1), numel(cols)])
    error('bench:runGedai:sizeMismatch', ...
        'GEDAI returned [%s] for [%d %d] input - samples were dropped or channels removed.', ...
        num2str(size(EEGclean.data)), size(data, 1), numel(cols));
end
Xhat(:, cols) = single(EEGclean.data);
end
