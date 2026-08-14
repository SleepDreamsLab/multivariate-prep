function failures = run_gedai_bids(BIDS, opts)
% RUN_GEDAI_BIDS  Run GEDAI artefact-rejection on pre-filtered BIDS EEG data.
%
%   run_gedai_bids(BIDS)
%   run_gedai_bids(BIDS, Name, Value, ...)
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
%                     run_filter_bids removed as bad.
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
end

if isempty(opts.inputpath), opts.inputpath = fullfile(BIDS.pth, 'derivatives', opts.derivfolder); end
if isempty(opts.savepath),  opts.savepath  = fullfile(BIDS.pth, 'derivatives', opts.derivfolder); end
if isempty(opts.figpath),   opts.figpath   = fullfile(BIDS.pth, 'derivatives', opts.derivfolder, 'figures'); end

KeepTime = struct();
if isempty(opts.runs), opts.runs = gedai.defaultRuns();
end

%%% Query EEG files from BIDS (used for entity extraction and subject iteration)
filesEEG = bids.query(BIDS, 'data', 'extension', '.vhdr', ...
    'task', opts.tasklabel, 'acq', opts.acqlabel);
if isempty(filesEEG); error('run_gedai_bids:noFiles', 'No matching EEG files found in BIDS layout.');
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
    
    %%% Try block
%     try

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
        error('run_gedai_bids:noScoring', 'No scoring file matched for %s.', fileID);
    end

    %%% Load sleep scoring
    fprintf('\nReading %s ...\n', scoringFile)
    scoringDigits = scoreloader(scoringFile);    

    %%% Import filtered EEG
    D = tic; fprintf('\nImporting filtered EEG ...\n')
    EEG = eeg_import(filtFile);
    KeepTime.EEGimport = toc(D);

    %%% Correct scoring length if needed
    nEpochs = floor(EEG.pnts / (opts.epochlength * EEG.srate));
    while numel(scoringDigits) > nEpochs; scoringDigits(end) = []; end

    %%% Drop non-EEG channels (in case filtered file still carries extras)
    EEG = pop_select(EEG, 'nochannel', intersect(1:EEG.nbchan, opts.noteegchannels));

    %%% Read SFP file (dome-solved channel locations)
    %%% Only for legacy inputs: run_filter_bids now writes .set, which already carries
    %%% chanlocs and urchanlocs. Re-reading them here would be wrong anyway, since it
    %%% assumes the file still holds the full montage in SFP order, and the bad
    %%% channels have already been dropped from it.
    if isfield(EEG, 'chanlocs') && ~isempty(EEG.chanlocs) && ...
            isfield(EEG.chanlocs, 'X') && ~isempty([EEG.chanlocs.X])
        fprintf('\nUsing channel locations stored in %s\n', opts.inputfileext)

    elseif ~isempty(opts.sfppath)
        sfpFile = gedai.matchSfpFile(opts.sfppath, p.entities.sub, p.entities.ses);
        fprintf('\nReading %s ...\n', sfpFile)
        chanlocs     = readlocs(sfpFile);
        chanlocs_reg = register_fiducials(chanlocs);
        EEG.chanlocs = chanlocs_reg(1:EEG.nbchan);

        % Urchanlocs
        EEG.urchanlocs = EEG.chanlocs;  
        for iCh = 1:numel(EEG.chanlocs)
            EEG.chanlocs(iCh).urchan = iCh;
        end

    elseif strcmp(BIDS.description.Name, {'ercp'})
        chanfile = fullfile(fileparts(eegFile), [fileID, '_channels.tsv']);
        elecfile = fullfile(fileparts(eegFile), ['sub-' p.entities.sub, '_ses-' p.entities.ses, '_electrodes.tsv']);
        [EEG, channelData, elecData] = bids_importchanlocs(EEG, chanfile, elecfile);

        % Urchanlocs
        EEG.urchanlocs = EEG.chanlocs;  
        for iCh = 1:numel(EEG.chanlocs)
            EEG.chanlocs(iCh).urchan = iCh;
        end
        
    else
        % continue
    end

    %%% Replace isolated N1 epochs at stage boundaries
    scoringDigits_NoN1 = gedai.killN1(scoringDigits);

    %%% Bad channels
    %%% Detected and removed by run_filter_bids before Zapline, where clean_channels
    %%% can still see the line noise its noise criterion is based on. The mask is only
    %%% needed here to select the matching rows of the leadfield.
    badchanFile = fullfile(opts.inputpath, subDir, [fileID '_desc-' opts.inputdesc '_badchans.mat']);
    if ~isfile(badchanFile)
        error('run_gedai_bids:noBadChannels', ...
            ['Bad channel file not found (%s). Bad channels are now detected and ' ...
             'removed in run_filter_bids, before Zapline - re-run it for this recording.'], badchanFile);
    end
    fprintf('Badchans → %s\n', badchanFile)
    removed_channels = getfield(load(badchanFile, 'removed_channels'), 'removed_channels');
    if numel(removed_channels) - nnz(removed_channels) ~= EEG.nbchan
        error('run_gedai_bids:badChannelMismatch', ...
            ['Bad channel mask expects %d channels left of %d, but the filtered file ' ...
             'has %d - mask and file are out of sync.'], ...
            numel(removed_channels) - nnz(removed_channels), numel(removed_channels), EEG.nbchan);
    end
    fprintf('%d channel(s) already removed as bad\n', nnz(removed_channels))

    %%% Leadfield covariance matrix
    lfCOV = gedai.loadrefcov(opts.leadfielddir, p, numel(removed_channels), removed_channels);

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
        gedaiFigDir  = fullfile(opts.figpath, ['desc-' opts.geddesc], fileID);
        gedaiDatFile = fullfile(gedaiRunDir, [fileID '_desc-' opts.geddesc '_eeg' opts.savefileext]);



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
            sidecarjson(KeepTime, ...
                fullfile(gedaiRunDir, [bvBase '.json']), ...
                struct('GEDAIParameters', rJson));
        end

        %%% Evaluation figures
        if ~exist(gedaiFigDir, 'dir'), mkdir(gedaiFigDir); end
        run.eval_clean(EEG, EEGgedai, scoringDigits_NoN1, ...
            'EpochLength', opts.epochlength, 'WelchWindow', 4, ...
            'EpochsToPlot', epochsToPlot, 'refresh', opts.refresh, ...
            'net', opts.net, 'SavePath', fullfile(gedaiFigDir, fileID));
        pause(1)
        close all;       
        
%         %%% Write savename marker
%         sidecarjson(KeepTime, ...
%             fullfile(gedaiFigDir, [savename '.json']), ...
%             struct('GEDAIParameters', r));

        %%% ICA residual
        if isfield(EEGgedai.etc, 'ic_classification')
            EEGgedai = ica.selectcomps(EEGgedai, 'ArtefactThreshold', 0.5, 'ManualQC', false);
            EEGgedai = pop_subcomp(EEGgedai, find(EEGgedai.reject.gcompreject), 0);
            run.eval_clean(EEG, EEGgedai, scoringDigits_NoN1, ...
                'EpochLength', opts.epochlength, 'WelchWindow', 4, ...
                'EpochsToPlot', epochsToPlot, 'refresh', opts.refresh, ...
                'SavePath', fullfile(gedaiFigDir, fileID));
            close all;
        end
    end

%     catch ME
%         fprintf('[ERROR] %s: %s\n', fileID, ME.message);
%         failures{end+1} = struct('fileID', fileID, 'message', ME.message, 'report', ME.getReport(), 'timestamp', datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')); %#ok<AGROW>
%     end
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

end % run_gedai_bids

% -------------------------------------------------------------------------
function mode = resolveStageMode(stages, dict)
% Resolve the cleaning mode for a group of stage digits from a stage-keyed dictionary.
% A stage group is cleaned in a single pass, so one strength has to cover all of it:
% stages within a group that disagree are an error, not something to silently pick from.
    if ~isConfigured(dict)
        error('run_gedai_bids:unconfiguredModeDict', ...
            'Mode dictionary is unconfigured - it must define a mode for every stage.');
    end

    missing = stages(~isKey(dict, stages));
    if ~isempty(missing)
        error('run_gedai_bids:missingStageMode', ...
            'No mode configured for stage(s) [%s] - add them to the mode dictionary.', ...
            num2str(missing(:)', '%d '));
    end

    modes = unique(string(dict(stages)));
    if numel(modes) > 1
        error('run_gedai_bids:inconsistentMode', ...
            ['Stages [%s] are cleaned together but map to different modes (%s) - ' ...
             'assign the same mode to all stages in a group.'], ...
            num2str(stages(:)', '%d '), strjoin(modes, ', '));
    end
    mode = char(modes);
end
