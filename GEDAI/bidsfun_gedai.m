function failures = bidsfun_gedai(BIDS, opts)
% BIDSFUN_GEDAI  Run GEDAI artefact-rejection on pre-filtered BIDS EEG data.
%
%   bidsfun_gedai(BIDS)
%   bidsfun_gedai(BIDS, Name, Value, ...)
%
%   Required
%   --------
%   BIDS              bids.layout object.
%
%   Input paths
%   -----------
%   filteredpath      Root of the filtered derivatives, organised as
%                     <sub>/<ses>/<fileID>_desc-<filtdesc>_eeg.set.
%                     Default: <BIDS root>/derivatives/prep-ged/A_filtered
%   filtdesc          desc label used when building the filtered filename.
%                     Default: 'filt'
%   inputfileext      Extension of the filtered input. Default: '.set' - the filtered
%                     file carries chanlocs/urchanlocs, which BrainVision cannot store
%                     and which are needed to interpolate the channels that
%                     bidsfun_hp_zap_cleanline removed as bad.
%   scoringpath       Directory containing sleep-scoring files (.json or .csv).
%                     Default: <BIDS root>/derivatives/scoring/scores/Manual_Checked
%   sfppath           Path passed to the SFP resolver.
%                     Default: <BIDS root>
%   leadfielddir      Root for Brainstorm leadfields (<sub>/<ses>/headmodel_surf_openmeeg.mat).
%                     Default: <BIDS root>/../Data_Analysis/Brainstorm_db/Leadfield_PM/data
%
%   Output paths
%   ------------
%   savepath          Root for GEDAI outputs.
%                     Default: <BIDS root>/derivatives/prep-ged/GEDAI
%   figpath           Root for all figures.
%                     Default: <BIDS root>/derivatives/prep-ged/figures
%   geddesc           desc label for GEDAI BrainVision output files.
%                     Default: 'filtGEDAI'
%   refresh           Force re-run even if a cache file exists. Default: false.
%
%   Parallel machines
%   -----------------
%   Several machines can be pointed at the same BIDS root and the same savepath at
%   once. Each recording is claimed with a lock file before any work starts, so a
%   second machine walking the same list skips whatever the first is busy with
%   instead of duplicating it. The claim is released when the recording finishes,
%   when it errors, and when the run is interrupted with Ctrl-C, so a failure never
%   parks a recording permanently. See claimFile.
%
%   uselocks          Claim each recording before processing it. Default: true.
%                     Set false for a single-machine run over a local savepath.
%   lockpath          Directory holding the claim files.
%                     Default: <savepath>/.locks
%   lockstalemin      Minutes after which a claim whose heartbeat stopped is taken
%                     over - the escape hatch for a machine that crashed or was
%                     rebooted mid-recording. Keep it well above the longest
%                     plausible single-recording runtime. Default: 360 (6 h).
%   lockheartbeatmin  Minutes between heartbeat writes on a held claim. Default: 5.
%
%   EEG
%   ---
%   tasklabel         BIDS task label(s) to query. Default: {'Sleep','sleep'}.
%   acqlabel    BIDS recording label to query. Default: '125Hz'.
%   noteegchannels    Channel indices to drop. Default: 257:300.
%   net               EEG net identifier passed to chans1020. Default: 'EGI256'.
%
%   Subject filter
%   --------------
%   subjectfilter     Cell array of subject ID strings; {} = all subjects.
%   sessionfilter     Cell array of session ID strings; {} = all sessions.
%
%   GEDAI
%   -----
%   runmode           'StageSpecific', 'StateWise', or 'WholeNight'. Default: 'StageSpecific'.
%   epochlength       Sleep-epoch length in seconds. Default: 30.
%   runs              Cell array of GEDAI run-config structs. Default: gedai.defaultRuns().
%   epochstoplot      Epoch indices for diagnostic figures. Default: auto.
%   prefix            Prefix prepended to the run savename. Default: ''.
%   dilaten           Epochs by which N2/N3/REM grow into neighbouring Wake/N1 epochs
%                     before killN1, so a boundary epoch holding real slow waves is
%                     cleaned as sleep rather than as wake. 0 disables. Default: 1.
%   dilatedirection   'forward' (default), 'backward' or 'both' - which side of a sleep
%                     run is claimed. See gedai.dilateStages.
%   pooltype          'Processes' (default) or 'Threads' for GEDAI's band loop. Threads
%                     share the broadcast data instead of copying it to every worker,
%                     which is where ~10% of pipeline runtime goes, but they also share
%                     one BLAS thread pool. Time it on one recording before adopting.
%
%   See dependancies.m.

arguments
    BIDS

    %--- Derivative folder ---
    opts.derivfolder      char = 'prep-ged'

    %--- Input paths ---
    opts.inputpath        char = ''
    opts.inputdesc        char = 'filt'
    opts.inputfileext     char = '.set'
    opts.scoringpath      char = fullfile(BIDS.pth, 'derivatives', 'scoring', 'scores', 'Manual_Checked')
    opts.sfppath          char = BIDS.pth
    opts.leadfielddir     char = fullfile(BIDS.pth, '..', 'Data_Analysis', 'Brainstorm_db', 'Leadfield_PM', 'data')

    %--- Output paths ---
    opts.savepath         char = ''
    opts.figpath          char = ''
    opts.geddesc          char = 'filt2ged'
    opts.refresh (1,1) logical = false
    opts.savefileext      char = '.set'

    %--- Parallel machines ---
    opts.uselocks (1,1) logical = true
    opts.lockpath         char  = ''
    opts.lockstalemin     (1,1) double {mustBePositive} = 360
    opts.lockheartbeatmin (1,1) double {mustBePositive} = 5

    %--- EEG ---
    opts.tasklabel                      = {'Sleep', 'sleep'}
    opts.acqlabel    char               = '125Hz'
    opts.noteegchannels    (1,:) double = 257:300
    opts.net               char         = 'EGI256'

    %--- Subject filter ---
    opts.subjectfilter     cell          = {}
    opts.sessionfilter     cell          = {}

    %--- GEDAI ---
    opts.runmode {mustBeMember(opts.runmode, {'WholeNight', 'StageSpecific', 'StateWise'})} = 'StageSpecific'
    opts.epochlength (1,1) double = 30
    opts.runs                     = []
    opts.epochstoplot             = []
    opts.prefix            char   = ''
    opts.pooltype {mustBeMember(opts.pooltype, {'Processes', 'Threads'})} = 'Processes'
    opts.evalplots (1,1) logical = true
    opts.dilaten (1,1) double {mustBeInteger, mustBeNonnegative} = 1
    opts.dilatedirection {mustBeMember(opts.dilatedirection, {'both','forward','backward'})} = 'forward'
end

fprintf('\n=== Running bidsfun_gedai ===\n');

if isempty(opts.inputpath), opts.inputpath = fullfile(BIDS.pth, 'derivatives', opts.derivfolder); end
if isempty(opts.savepath),  opts.savepath  = fullfile(BIDS.pth, 'derivatives', opts.derivfolder); end
if isempty(opts.figpath),   opts.figpath   = fullfile(BIDS.pth, 'derivatives', opts.derivfolder, 'figures'); end
if isempty(opts.lockpath),  opts.lockpath  = fullfile(opts.savepath, '.locks'); end

KeepTime = struct();
if isempty(opts.runs), opts.runs = gedai.defaultRuns();
end

%%% Query EEG files from BIDS (used for entity extraction and subject iteration)
filesEEG = bids.query(BIDS, 'data', 'extension', '.vhdr', ...
    'task', opts.tasklabel, 'acq', opts.acqlabel);
if isempty(filesEEG); error('bidsfun_gedai:noFiles', 'No matching EEG files found in BIDS layout.');
end

%%% Scoring files
if ~isempty(opts.scoringpath); scoringfiles = gedai.collectScoringFiles(opts.scoringpath);
end

%%% Loop over EEG files
failures = {};
for ifile = 1:numel(filesEEG)
    eegFile = filesEEG{ifile};
    p       = bids.internal.parse_filename(eegFile);
    fileID  = strjoin(cellfun(@(k) [k '-' p.entities.(k)], fieldnames(p.entities), 'uni', 0), '_');
    subDir  = fullfile(['sub-' p.entities.sub], ['ses-' p.entities.ses]);

    %%% Subject filter
    if ~isempty(opts.subjectfilter) && ~contains(fileID, opts.subjectfilter)
        continue
    end

    %%% Session filter
    if ~isempty(opts.sessionfilter) && ~contains(fileID, opts.sessionfilter)
        continue
    end
    fprintf('\n=== %s ===\n', fileID)

    %%% Skip this file entirely if every configured run's output already exists and
    %%% refresh is not requested.
    if ~opts.refresh && ~runsPending(opts, subDir, fileID)
        fprintf('[skip] all run output(s) exist for %s\n', fileID)
        continue
    end

    %%% Claim this recording, so a second machine walking the same list moves on to the
    %%% next one instead of redoing this. The lock file is created with an atomic
    %%% create-if-absent, so two machines reaching this line together cannot both win.
    %%% lockGuard holds the claim: it releases on success, on the error caught below,
    %%% and on Ctrl-C or any error that unwinds out of this function, so a recording is
    %%% never left claimed by a run that is no longer working on it. Hence the explicit
    %%% clear at the end of the iteration - the claim must not outlive its recording.
    lockGuard = [];  %#ok<NASGU> release any claim still held from the previous iteration
    if opts.uselocks
        [lockGuard, acquired, holder] = claimFile( ...
            fullfile(opts.lockpath, [fileID '.lock']), ...
            'stalemin', opts.lockstalemin, 'heartbeatmin', opts.lockheartbeatmin); %#ok<ASGLU>
        if ~acquired
            fprintf('[skip] %s: claimed by %s (pid %d) since %s\n', ...
                fileID, holder.host, holder.pid, holder.started)
            continue
        end
        %%% Re-check now that the claim is ours: the other machine may have finished
        %%% this recording and dropped its lock between our skip check above and here.
        if ~opts.refresh && ~runsPending(opts, subDir, fileID)
            fprintf('[skip] all run output(s) exist for %s\n', fileID)
            continue
        end
    end

    %%% Try block
    try

        %%% Resolve filtered input file
        filtFile = fullfile(opts.inputpath, subDir, [fileID '_desc-' opts.inputdesc '_eeg' opts.inputfileext]);
        if ~isfile(filtFile)
            fprintf('[skip] %s: filtered file not found (%s)\n', fileID, filtFile)
            continue
        end
        fprintf('Input  → %s\n', filtFile)
        fprintf('Output → %s\n', opts.savepath)

        %%% Find matching scoring file
        if isempty(opts.scoringpath)
            scoringpath = fullfile(BIDS.pth, subDir);
            scoringfiles = gedai.collectScoringFiles(scoringpath);
        end
        scoringFile = gedai.matchScoringFile(p.entities, scoringfiles);
        fprintf('Scoring → %s\n', scoringFile)
        if isempty(scoringFile)
            error('bidsfun_gedai:noScoring', 'No scoring file matched for %s.', fileID);
        end

        %%% Load sleep scoring
        fprintf('\nReading %s ...\n', scoringFile)
        scoringDigits = scoreloader(scoringFile);

        %%% Import filtered EEG
        D = tic; fprintf('\nImporting filtered EEG ...\n')
        EEG = fast_eeg_import(filtFile);
        KeepTime.EEGimport = toc(D);

        %%% Correct scoring length if needed
        nEpochs = floor(EEG.pnts / (opts.epochlength * EEG.srate));
        while numel(scoringDigits) > nEpochs; scoringDigits(end) = []; end

        %%% Drop non-EEG channels (in case filtered file still carries extras)
        elabels = arrayfun(@(x) ['E' num2str(x)], opts.noteegchannels, 'uni', 0);
        EEG = pop_select(EEG, 'nochannel', intersect({EEG.chanlocs.labels}, elabels));

        %%% Channel locations
        %%% Only for legacy inputs: bidsfun_hp_zap_cleanline now writes .set, which already
        %%% carries chanlocs and urchanlocs. gedai.assignChanlocs leaves those alone and only
        %%% falls back to the SFP/ercp sources when they are missing.
        EEG = gedai.assignChanlocs(EEG, BIDS, opts.sfppath, eegFile, p, fileID);

        %%% Grow the sleep stages into the following epoch, then resolve N1.
        %%% Scoring labels a whole 30-s epoch, so the epoch after a sleep run can still hold
        %%% genuine slow waves; cleaning it as Wake risks removing them. Done before killN1
        %%% so that an N1 epoch touching sleep is claimed by that sleep stage first, which
        %%% also shifts killN1's split of any longer N1 run one epoch towards sleep.
        %%% Cleaning-time only - scoringDigits itself is untouched.
        [scoringDilated, nDilated] = gedai.dilateStages(scoringDigits, ...
            'n', opts.dilaten, 'direction', opts.dilatedirection);
        scoringDigits_NoN1 = gedai.killN1(scoringDilated);

        %%% Bad channels
        %%% bidsfun_hp_zap_cleanline already removed them before saving the .set, and every
        %%% surviving channel's urchan points back into EEG.urchanlocs (both are set before
        %%% any channel is ever dropped, in bidsfun_hp_zap_cleanline). So which channels are missing -
        %%% and therefore which rows of the leadfield to use - is already in the file; no
        %%% separate mask file needs loading.
        removed_channels = true(numel(EEG.urchanlocs), 1);
        removed_channels([EEG.chanlocs.urchan]) = false;
        if nnz(~removed_channels) ~= EEG.nbchan
            error('bidsfun_gedai:badChannelMismatch', ...
                ['urchan indices on EEG.chanlocs do not match EEG.urchanlocs (%d of %d ' ...
                'channels resolved) - chanlocs and urchanlocs are out of sync.'], ...
                nnz(~removed_channels), EEG.nbchan);
        end
        fprintf('%d channel(s) were removed as bad previously\n', nnz(removed_channels))

        %%% Leadfield covariance matrix
        fprintf('removing the same %d bad channels, as well as non-EEG channels from leadfield matrix ...\n', ...
            nnz(removed_channels))
        lfCOV = gedai.loadrefcov(opts.leadfielddir, p, numel(removed_channels), removed_channels, opts.noteegchannels);

        %%% Average re-reference
        EEG.data = EEG.data - sum(EEG.data, 1) / (size(EEG.data, 1) + 1);

        %%% Epochs to plot
        epochsToPlot = gedai.resolveEpochsToPlot(opts.epochstoplot, scoringDigits);

        %%% GEDAI runs
        for iRun = 1:numel(opts.runs)
            r = opts.runs{iRun};

            switch opts.runmode
                case 'StageSpecific'
                    stageLogic      = {[-2], [-3], [0], [1]};
                    refCOV_perStage = {lfCOV, lfCOV, lfCOV, lfCOV};
                case 'StateWise'
                    stageLogic      = {[-2, -3], [0, 1]};
                    refCOV_perStage = {lfCOV, lfCOV};
                case 'WholeNight'
                    % N1 (-1) is deliberately excluded: killN1 (called earlier) reassigns every
                    % N1 epoch to a neighbouring stage, so -1 never appears in scoringDigits_NoN1.
                    % Using the -3:1 range here would still ask the mode dict for a -1 entry.
                    stageLogic      = {[-3, -2, 0, 1]};
                    refCOV_perStage = {lfCOV};
            end

            %%% Build per-stage GEDAIMode / GEDAIModeBB from stage-keyed dicts
            GEDAIMode_perStage   = cellfun(@(s) resolveStageMode(s, r.GEDAIMode_dict),   stageLogic, 'uni', 0);
            GEDAIModeBB_perStage = cellfun(@(s) resolveStageMode(s, r.GEDAIModeBB_dict), stageLogic, 'uni', 0);

            %%% Run name: built from the resolved modes, and from the full stage list
            %%% before recording-specific stages are dropped, so it describes the
            %%% configuration and stays comparable across recordings
            savename = [opts.prefix opts.runmode '_' ...
                gedai.buildSaveName(r, EEG.srate, GEDAIMode_perStage, GEDAIModeBB_perStage)];
            
            geddesc = opts.geddesc;
            if isempty(opts.geddesc)
                geddesc = savename;
            end

            %%% Drop stage groups absent from scoring
            presentStages      = unique(scoringDigits_NoN1);
            keep               = cellfun(@(s) any(ismember(s, presentStages)), stageLogic);
            stageLogic           = stageLogic(keep);
            refCOV_perStage      = refCOV_perStage(keep);
            GEDAIMode_perStage   = GEDAIMode_perStage(keep);
            GEDAIModeBB_perStage = GEDAIModeBB_perStage(keep);

            fprintf('Run %d/%d: %s\n', iRun, numel(opts.runs), savename)

            %%% Paths for this run
            gedaiRunDir  = fullfile(opts.savepath, subDir);
            gedaiFigDir  = fullfile(opts.figpath, ['desc-' geddesc], fileID);
            gedaiDatFile = fullfile(gedaiRunDir, [fileID '_desc-' geddesc '_eeg' opts.savefileext]);

            %%% Skip run if output already exists and refresh not requested
            if ~opts.refresh && isfile(gedaiDatFile)
                fprintf('[skip] output exists: %s\n', gedaiDatFile)
                continue
            end

            %%% Run GEDAI (cached as BrainVision .dat)
            clear EEGgedai
            isNewRun = opts.refresh || ~isfile(gedaiDatFile);
            D = tic;
            EEGgedai = smartcache( ...
                @() run.GEDAI_StageSpecific(EEG, scoringDigits_NoN1, ...
                stageLogic, KeepTime, ...
                'EpochLength',                opts.epochlength, ...
                'GEDAIMode',                  GEDAIMode_perStage, ...
                'GEDAIModeBB',                GEDAIModeBB_perStage, ...
                'GEDAIEpochSize',             r.GEDAIEpochSize, ...
                'GEDAILowCutOffFreq',         r.GEDAILowCutOffFreq, ...
                'BBEpochSize',                r.GEDAIBroadbandEpochSize, ...
                'BroadbandOnly',              r.broadbandOnly, ...
                'GEDAIEnovaChannelThreshold', r.GEDAIEnovaChannelThreshold, ...
                'PercentileThreshold',        r.percentileThreshold, ...
                'BBMinThreshold',             r.BBMinThreshold, ...
                'ComputeSENSAI',              r.computeSENSAI, ...
                'ICAtype',                    r.ICAtype, ...
                'PoolType',                   opts.pooltype, ...
                'RefCOV',                     refCOV_perStage), ...
                gedaiDatFile, opts.refresh, {'EEGgedai', '', '', ''});
            KeepTime.GEDAI = toc(D);
            fprintf('GEDAI took %.2f min\n', KeepTime.GEDAI / 60)

            %%% JSON sidecar: GEDAI parameters + timing (only written on fresh runs)
            if isNewRun
                [~, bvBase] = fileparts(gedaiDatFile);
                rJson = rmfield(r, {'GEDAIMode_dict', 'GEDAIModeBB_dict'});
                rJson.GEDAIMode_resolved   = GEDAIMode_perStage;
                rJson.GEDAIModeBB_resolved = GEDAIModeBB_perStage;
                rJson.StageDilation = struct('n', opts.dilaten, ...
                    'direction', opts.dilatedirection, 'nEpochsRelabelled', nDilated);
                sidecarjson(KeepTime, ...
                    fullfile(gedaiRunDir, [bvBase '.json']), ...
                    struct('GEDAIParameters', rJson));
            end

            %%% Evaluation figures. Optional: bidsfun_evalfigs can produce the same
            %%% comparison as a separate stage, which is the better route when the figures
            %%% are being re-made without re-running GEDAI.
            if opts.evalplots
                if ~exist(gedaiFigDir, 'dir'), mkdir(gedaiFigDir); end
                run.eval_clean(EEG, EEGgedai, scoringDigits_NoN1, ...
                    'EpochLength', opts.epochlength, 'WelchWindow', 4, ...
                    'EpochsToPlot', epochsToPlot, 'refresh', opts.refresh, ...
                    'net', opts.net, 'SavePath', fullfile(gedaiFigDir, fileID));
                pause(1)
                close all;
                pause(4)
            end

            %         %%% Write savename marker
            %         sidecarjson(KeepTime, ...
            %             fullfile(gedaiFigDir, [savename '.json']), ...
            %             struct('GEDAIParameters', r));

            %%% ICA residual
            if opts.evalplots && isfield(EEGgedai.etc, 'ic_classification')
                EEGgedai = ica.selectcomps(EEGgedai, 'ArtefactThreshold', 0.5, 'ManualQC', false);
                EEGgedai = pop_subcomp(EEGgedai, find(EEGgedai.reject.gcompreject), 0);
                run.eval_clean(EEG, EEGgedai, scoringDigits_NoN1, ...
                    'EpochLength', opts.epochlength, 'WelchWindow', 4, ...
                    'EpochsToPlot', epochsToPlot, 'refresh', opts.refresh, ...
                    'SavePath', fullfile(gedaiFigDir, fileID));
                close all;
            end
        end

    catch ME
        fprintf('[ERROR] %s: %s\n', fileID, ME.message);
        failures{end+1} = struct('fileID', fileID, 'message', ME.message, 'report', ME.getReport(), 'timestamp', datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')); %#ok<AGROW>
    end

    %%% Drop the claim, whether the recording succeeded or failed. A failed recording
    %%% has to become available again - it is exactly the one another machine (or a
    %%% later run of this one) should be free to retry.
    clear lockGuard
end

%%% Failure summary
if ~isempty(failures)
    fprintf('\n=== %d file(s) failed ===\n', numel(failures));
    for k = 1:numel(failures)
        fprintf('  %s: %s\n', failures{k}.fileID, failures{k}.message);
    end
    if ~exist(opts.savepath, 'dir'), mkdir(opts.savepath); end
    fid = fopen(fullfile(opts.savepath, 'failed_files_gedai.json'), 'w');
    fprintf(fid, '%s', jsonencode([failures{:}], 'PrettyPrint', true));
    fclose(fid);
end

end % bidsfun_gedai

% -------------------------------------------------------------------------
function tf = runsPending(opts, subDir, fileID)
% True when at least one configured run has no output yet. Checked before the scoring
% load and the (slow) EEG import, since none of that is needed just to decide there is
% nothing to do. The per-run check further down stays as well, to still resume
% correctly if only some runs' outputs are missing.
tf = false;
for iRun = 1:numel(opts.runs)
    gedaiDatFile = fullfile(opts.savepath, subDir, ...
        [fileID '_desc-' opts.geddesc '_eeg' opts.savefileext]);
    if ~isfile(gedaiDatFile)
        tf = true;
        return
    end
end
end

% -------------------------------------------------------------------------
function mode = resolveStageMode(stages, dict)
% Resolve the cleaning mode for a group of stage digits from a stage-keyed dictionary.
% A stage group is cleaned in a single pass, so one strength has to cover all of it:
% stages within a group that disagree are an error, not something to silently pick from.
if ~isConfigured(dict)
    error('bidsfun_gedai:unconfiguredModeDict', ...
        'Mode dictionary is unconfigured - it must define a mode for every stage.');
end

missing = stages(~isKey(dict, stages));
if ~isempty(missing)
    error('bidsfun_gedai:missingStageMode', ...
        'No mode configured for stage(s) [%s] - add them to the mode dictionary.', ...
        num2str(missing(:)', '%d '));
end

modes = unique(string(dict(stages)));
if numel(modes) > 1
    error('bidsfun_gedai:inconsistentMode', ...
        ['Stages [%s] are cleaned together but map to different modes (%s) - ' ...
        'assign the same mode to all stages in a group.'], ...
        num2str(stages(:)', '%d '), strjoin(modes, ', '));
end
mode = char(modes);
end
