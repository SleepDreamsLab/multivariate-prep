function [EEGclean, EEGstage, epochIdx, KeepTime] = GEDAI_PerStage(EEG, Scoring, SleepStages, KeepTime, opts)
% GEDAI_PERSTAGE  Run GEDAI separately on each sleep-stage group and merge.
%
%   Inputs
%   ------
%   EEG         : EEGLAB struct (continuous or pre-epoched).
%   Scoring     : Numeric vector of stage digits, one per epoch.
%   SleepStages : Cell array — each cell lists the stage digit(s) for one
%                 GEDAI run, e.g. {[-2 -3], [0], [1]}.  Default: {[0]}.
%   KeepTime    : Timing struct; pass [] to start fresh.
%
%   Name-value options
%   ------------------
%   EpochLength        : Epoch length in seconds (default 30).
%   GEDAIMode          : Cell array of mode strings, one per SleepStages group.
%                        If shorter than nGroups, the last entry is recycled.
%                        Default: {'auto-'}.
%   GEDAIEpochSize     : Epoch size in cycles (default 12).
%   GEDAILowCutOffFreq : Low-cut frequency Hz (default 0.1).
%   GEDAIMethod        : Reference string used when RefCOV is absent
%                        (default 'interpolated').
%   RefCOV             : Reference covariance matrix — three forms accepted:
%                          cell array  — one entry per SleepStages group;
%                                        [] in a cell → GEDAIMethod for
%                                        that group.
%                          matrix      — same matrix used for all groups.
%                          [] / {}     — all groups use GEDAIMethod.
%   MovAvgSize         : GEDAI moving-average size (default 60).
%   BBEpochSize        : GEDAI broadband epoch size (default 12).
%
%   Outputs
%   -------
%   EEGclean  : GEDAI-cleaned continuous EEG, all groups merged and
%               re-sorted into the original epoch order.
%   EEGstage  : Uncleaned EEG for the same epochs.
%   epochIdx  : Epoch indices (into Scoring) for all processed epochs.
%   KeepTime  : Input struct extended with fields Epoching and GEDAI (s).

arguments
    EEG
    Scoring
    SleepStages             = {[0]}
    KeepTime                = []
    opts.EpochLength (1,1)  = 30
    opts.GEDAIMode          = {'auto-'}
    opts.GEDAIModeBB        = {'auto-'}
    opts.GEDAIEpochSize     = 12
    opts.GEDAILowCutOffFreq = 0.1
    opts.GEDAIMethod        = 'interpolated'
    opts.RefCOV             = {}    % cell | matrix | [] — see above
    opts.MovAvgSize         = []
    opts.BBEpochSize        = 12
    opts.BroadbandOnly (1,1) logical = false
    opts.PercentileThreshold        = []   % [] → auto-derived inside GEDAI
    opts.BBMinThreshold             = -2   % broadband minimum threshold
    opts.ComputeSENSAI (1,1) logical = true
    opts.ICAtype                     = 'none'
    opts.GEDAIEnovaChannelThreshold = Inf;
    opts.PoolType {mustBeMember(opts.PoolType, {'Processes', 'Threads'})} = 'Processes'
end

if ~iscell(SleepStages)
    SleepStages = {SleepStages};
end

Scoring = Scoring(:)';

% ── 1. Epoch continuous data ──────────────────────────────────────────────
if EEG.trials == 1
    D = tic;
    % Strip boundary events — they cause eeg_regepochs to silently drop
    % overlapping epochs, creating a mismatch with Scoring.
    if ~isempty(EEG.event)
        EEG.event(strcmpi({EEG.event.type}, 'boundary')) = [];
        EEG.event(contains({EEG.event.type}, 'Epoch')) = [];
        EEG = eeg_checkset(EEG);
    end
    EEG = eeg_regepochs(EEG, 'recurrence', opts.EpochLength, ...
        'limits', [0 opts.EpochLength], ...
        'eventtype', sprintf('Epoch%ds', opts.EpochLength));
    KeepTime.Epoching = toc(D);
end

if numel(Scoring) ~= EEG.trials
    error('gedai.GEDAI_PerStage: Scoring length (%d) ≠ number of epochs (%d).', ...
        numel(Scoring), EEG.trials);
end

% ── 2. Loop over stage groups ─────────────────────────────────────────────
nGroups       = numel(SleepStages);
EEGclean_list = cell(nGroups, 1);
EEGstage_list = cell(nGroups, 1);
epochIdx_list = cell(nGroups, 1);
gedaiMetrics  = struct([]);
gedaiTime     = 0;

% % EGI-256 → 10-20 mapping (used for label renaming after GEDAI)
% egi256_to_1020 = struct( ...
%     'Fp1',  37, 'Fpz',  26, 'Fp2',  18, ...
%     'F7',   47, 'F3',   36, 'Fz',   21, 'F4', 224, 'F8',   2, ...
%     'T7',   69, 'C3',   59, 'C4',  183, 'T8', 202, ...
%     'P7',   96, 'P3',   87, 'Pz',  101, 'P4', 153, 'P8', 170, ...
%     'O1',  116, 'O2',  150);
% egi256indices = struct2array(egi256_to_1020);
% labels1020    = fieldnames(egi256_to_1020);

for iGroup = 1:nGroups
    stages = SleepStages{iGroup};
    idx    = find(ismember(Scoring, stages));

    if isempty(idx)
        error('gedai.GEDAI_PerStage: no epochs found for stage(s) [%s].', ...
            num2str(stages(:)', '%d '));
    end

    fprintf('gedai.GEDAI_PerStage: group %d/%d – %d epochs, stage(s) [%s] (%.1f min).\n', ...
        iGroup, nGroups, numel(idx), num2str(stages(:)', '%d '), ...
        numel(idx) * opts.EpochLength / 60);

    epochIdx_list{iGroup} = idx;

    EEGstageGroup = pop_select(EEG, 'trial', idx);
    EEGstageGroup = eeg_epoch2continuous(EEGstageGroup);

    % Resolve reference matrix for this group
    if iscell(opts.RefCOV) && iGroup <= numel(opts.RefCOV) && ~isempty(opts.RefCOV{iGroup})
        gedaiRefMatrix = opts.RefCOV{iGroup};       % per-group matrix
    elseif isnumeric(opts.RefCOV) && ~isempty(opts.RefCOV)
        gedaiRefMatrix = opts.RefCOV;               % single matrix → all groups
    else
        gedaiRefMatrix = opts.GEDAIMethod;          % string fallback
    end

    % Resolve GEDAIMode / GEDAIModeBB for this group (recycle last entry if cell is shorter)
    gedaiMode   = opts.GEDAIMode{min(iGroup,   numel(opts.GEDAIMode))};
    gedaiModeBB = opts.GEDAIModeBB{min(iGroup, numel(opts.GEDAIModeBB))};
    fprintf('Running GEDAI with mode=%s, modeBB=%s\n', gedaiMode, gedaiModeBB)



    % ── Worker count from available RAM. Peak per-worker allocation in
    %    GEDAI_per_band is roughly the band's data copy plus working
    %    copies; BYTES_PER_SAMPLE_CH bundles the double (8 B) with that
    %    multiplier.
    %
    %    Recalibrated 2026-08-13 after removing GEDAI_per_band's per-band
    %    COV/COV_2/Eval/Eval_2 arrays and rank-truncating the GEVD (see
    %    GEDAI-svennonito commit 8c39921). Measured directly: one worker's
    %    full peak footprint (broadcast copy of unfiltered_data + its
    %    MODWT band extraction + its GEDAI_per_band call), N=256, at three
    %    matched-scale bands (T=48/192/384 samples), P swept 450k-1.8M
    %    samples to confirm linearity in N*P before trusting the ratio:
    %       old code : 46.4 x (N*P*8B)   [band 1, T=48 — matches the report's
    %                   dominant-band claim; note this is already ~2x the
    %                   24x hardcoded above, so worker counts before this
    %                   change may have been under-provisioned even pre-fix]
    %       new code : 15.2-15.3 x (N*P*8B), consistent across all three
    %                   bands tested (T=48, 192, 384) to within 1%. The
    %                   per-band GEVD arrays that used to set the ceiling
    %                   are now small enough that the fixed, band-independent
    %                   costs (broadcast copy, MODWT scratch, clean_EEG's
    %                   output buffers) dominate instead — hence no strong
    %                   per-band variation left to calibrate against.
    %    16x carried ~5% headroom over that measured ceiling; the pre-existing
    %    SAFETY factor below covers the rest. Re-measure if channel count,
    %    MODWT extraction (e.g. porting stateful_modwt_single_band), or
    %    clean_EEG's output buffering changes again.
    %
    %    Halved to 8x on switching GEDAI to parallel = 'blocks'. The 15-16x
    %    figure was dominated by each worker holding a whole band; under block
    %    parallelism a worker holds one time block instead, so the per-worker
    %    peak is a fraction of the band. 8x is a deliberately conservative
    %    stand-in for that - it has NOT been measured the way the 15.2x was, so
    %    if a night starts swapping this is the first number to re-measure.
    BYTES_PER_SAMPLE_CH = 8 * 1;    % 8 B double x ~8 working copies (blocks)
    SAFETY              = 0.99;      % leave headroom for client + OS cache
    %%% GEDAI is called with parallel = 'blocks', so the band loop runs serially and the
    %%% parfor is over time blocks within each band. That removes the old ceiling: there
    %%% are far more blocks than the 11 wavelet bands, and they are equal-sized, so extra
    %%% workers keep earning instead of idling. (The band-count cap that used to live
    %%% here applied to parallel = true, where the loop had one task per band.)
    %%% The ceiling is now whatever the local cluster profile allows. Inf does not work:
    %%% the memory formula below can exceed the profile's NumWorkers, and parpool errors
    %%% when asked for more workers than the cluster has. Raise it in the Parallel
    %%% preferences ("Processes" profile) if you want more than this reports.
    localCluster = parcluster('Processes');
    MAX_WORKERS  = localCluster.NumWorkers;

    LOAD_NOW = EEGstageGroup.pnts * EEGstageGroup.nbchan;
    perWorkerBytes = LOAD_NOW * BYTES_PER_SAMPLE_CH;

    [~, sysMem] = memory;                       % Windows only
    freeBytes   = sysMem.PhysicalMemory.Available;

    nWorkers = floor(SAFETY * freeBytes / perWorkerBytes);
    nWorkers = max(1, min(MAX_WORKERS, nWorkers));

    %%% Pool type. GEDAI's band loop broadcasts unfiltered_data to every worker; on a
    %%% process pool that is a full serialized copy each (~7 GB for 4 h of N2 at 256 ch),
    %%% and the profile attributes ~10% of pipeline runtime to that dispatch. A thread
    %%% pool shares the array instead of copying it. The trade is that thread workers
    %%% share one BLAS thread pool, so per-band svd/eig calls contend where separate
    %%% processes would not - hence an option to be timed, not a new default.
    wantThreads = strcmpi(opts.PoolType, 'Threads');
    if wantThreads
        %%% A thread pool cannot exceed MATLAB's computational threads, and the budget
        %%% above assumes a per-worker copy that threads do not make, so it
        %%% under-provisions here rather than over.
        nWorkers = max(1, min(nWorkers, maxNumCompThreads));
    end

    % Build target pool (also rebuild when the pool that is up is of the wrong type)
    p             = gcp('nocreate');
    poolIsThreads = ~isempty(p) && contains(class(p), 'ThreadPool');
    if isempty(p) || poolIsThreads ~= wantThreads || abs(p.NumWorkers - nWorkers) > 0
        delete(gcp('nocreate'));
        if wantThreads
            parpool('Threads', nWorkers);
        else
            c = localCluster;       % same object the worker cap was read from
            oldNT = c.NumThreads;
            c.NumThreads = 1;
            parpool(c, nWorkers);
            c.NumThreads = oldNT;   % restore profile
        end
    end

    fprintf(['gedai.GEDAI_PerStage: load %.2e, %.1f GB/worker, ' ...
             '%.1f GB free -> %d workers.\n'], ...
        LOAD_NOW, perWorkerBytes/2^30, freeBytes/2^30, nWorkers);





    % Run GEDAI
    D = tic;
    fprintf('gedai.GEDAI_PerStage: running GEDAI for stage(s) [%s] ...\n', ...
        num2str(stages(:)', '%d '));
    [EEGcleanGroup, ~, SENSAI_score, SENSAI_score_per_band, ...
        artifact_threshold_per_band, mean_ENOVA, ENOVA_per_epoch] = ...
        GEDAI(EEGstageGroup, gedaiMode, gedaiModeBB, opts.GEDAIEpochSize, ...
              opts.GEDAILowCutOffFreq, gedaiRefMatrix, 'blocks', 0, [], opts.GEDAIEnovaChannelThreshold, [], ...
              opts.MovAvgSize, opts.BBEpochSize, opts.BroadbandOnly, opts.PercentileThreshold, opts.BBMinThreshold, opts.ComputeSENSAI);
    gedaiTime = gedaiTime + toc(D);

    % Store per-group metrics
    EEGcleanGroup.etc.GEDAI.sleepStage      = stages; 
    EEGcleanGroup.etc.GEDAI.epochIdx        = idx; 
    EEGcleanGroup.etc.GEDAI.time_spent_min  = round(toc(D) / 60, 2); 

%     % Store per-group metrics
%     m.sleepStage                  = stages;
%     m.epochIdx                    = idx;
% %     m.SENSAI_score                = SENSAI_score;
% %     m.SENSAI_score_per_band       = SENSAI_score_per_band;
% %     m.artifact_threshold_per_band = artifact_threshold_per_band;
% %     m.mean_ENOVA                  = mean_ENOVA;
% %     m.ENOVA_per_epoch             = ENOVA_per_epoch;
%     m.time_spent_min              = round(toc(D) / 60, 2);
%     if isempty(gedaiMetrics); gedaiMetrics = m;
%     else;                     gedaiMetrics(end+1) = m; end %#ok<AGROW>

%     % Rename 10-20 channel labels in-place (does not reduce channel count)
%     for iCh = 1:numel(egi256indices)
%         EEGcleanGroup.chanlocs(egi256indices(iCh)).labels = labels1020{iCh};
%         EEGstageGroup.chanlocs(egi256indices(iCh)).labels = labels1020{iCh};
%     end

    EEGclean_list{iGroup} = EEGcleanGroup;
    EEGstage_list{iGroup} = EEGstageGroup;
end

KeepTime.GEDAI = gedaiTime;

% ── 3. Merge all groups ───────────────────────────────────────────────────
epochIdx = cat(2, epochIdx_list{:});

EEGclean = EEGclean_list{1};
EEGstage = EEGstage_list{1};
G        = EEGclean_list{1}.etc.GEDAI;
for iGroup = 2:nGroups
    EEGclean = pop_mergeset(EEGclean, EEGclean_list{iGroup});
    EEGstage = pop_mergeset(EEGstage, EEGstage_list{iGroup});
    G(end+1) = EEGclean_list{iGroup}.etc.GEDAI;
end

%%% GEDAI Info
EEGclean.etc.GEDAI = G;


% ── 4. Restore original epoch order ──────────────────────────────────────
[epochIdx, sortOrder] = sort(epochIdx);
EEGclean = utils.resequence_epochs(EEGclean, sortOrder);
EEGstage = utils.resequence_epochs(EEGstage, sortOrder);

% EEGclean.etc.gedaiMeta = gedaiMetrics;

% ── 5. Optional ICA ───────────────────────────────────────────────────────
if ~strcmpi(opts.ICAtype, 'none')
    D = tic;
    fprintf('gedai.GEDAI_PerStage: running ICA (%s) once on whole night ...\n', opts.ICAtype);
    EEGclean = pop_runica(EEGclean, 'icatype', opts.ICAtype);
    EEGclean = iclabel(EEGclean);
    KeepTime.ICA = toc(D);
else
    KeepTime.ICA = 0;
end

EEGclean.data = single(EEGclean.data);
EEGstage.data = single(EEGstage.data);

end
