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
%                     <sub>/<ses>/<fileID>_desc-<filtdesc>_eeg.dat.
%                     Default: <BIDS root>/derivatives/prep-ged/A_filtered
%   filtdesc          desc label used when building the filtered filename.
%                     Default: 'filt'
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
%   noteegchannels    Channel indices to drop. Default: 257:264.
%   net               EEG net identifier passed to chans1020. Default: 'EGI256'.
%
%   Subject filter
%   --------------
%   subjectfilter     Cell array of subject ID strings; {} = all subjects.
%
%   GEDAI
%   -----
%   runmode           'StageSpecific', 'StateWise', or 'WholeNight'. Default: 'StageSpecific'.
%   epochlength       Sleep-epoch length in seconds. Default: 30.
%   runs              Cell array of GEDAI run-config structs. Default: gedai.defaultRuns().
%   epochstoplot      Epoch indices for diagnostic figures. Default: auto.
%   prefix            Prefix prepended to the run savename. Default: ''.
%
%   See dependancies.m.

arguments
    BIDS

    %--- Derivative folder ---
    opts.derivfolder      char = 'prep-ged'

    %--- Input paths ---
    opts.inputpath        char = ''
    opts.inputdesc        char = 'filt'
    opts.inputfileext     char = '.vhdr'    
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

    %--- GEDAI ---
    opts.runmode {mustBeMember(opts.runmode, {'WholeNight', 'StageSpecific', 'StateWise'})} = 'StageSpecific'
    opts.epochlength (1,1) double = 30
    opts.runs                     = []
    opts.epochstoplot             = []
    opts.prefix            char   = ''
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
    fprintf('\n=== %s ===\n', fileID)
    
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
    if ~isempty(opts.sfppath)
        sfpFile = gedai.matchSfpFile(opts.sfppath, p.entities.sub, p.entities.ses);
        fprintf('\nReading %s ...\n', sfpFile)
        chanlocs     = readlocs(sfpFile);
        chanlocs_reg = register_fiducials(chanlocs);
        EEG.chanlocs = chanlocs_reg(1:EEG.nbchan);
        EEG.urchanlocs = EEG.chanlocs;    
    elseif strcmp(BIDS.description.Name, {'ercp'})
        chanfile = fullfile(fileparts(rawFile), [fileID, '_channels.tsv']);
        elecfile = fullfile(fileparts(rawFile), ['sub-' p.entities.sub, '_ses-' p.entities.ses, '_electrodes.tsv']);
        [EEG, channelData, elecData] = bids_importchanlocs(EEG, chanfile, elecfile);
    else
        continue
        EEG.urchanlocs = EEG.chanlocs;            
    end

    %%% Average re-reference
    EEG.data = EEG.data - sum(EEG.data, 1) / (size(EEG.data, 1) + 1);

    %%% Replace isolated N1 epochs at stage boundaries
    scoringDigits_NoN1 = gedai.killN1(scoringDigits);

    %%% Assign 10-20 labels
    EEG = chans1020(EEG, false, 'net', opts.net);

    %%% Bad channel detection
    D = tic;
    [removed_channels, corr, znoise] = smartcache( ...
        @() clean_channels(EEG, 0.7, 4, [], 0.5, 25), ...
        fullfile(opts.savepath, subDir, [fileID '_badchans.mat']), ...
        false, {'', 'removed_channels', 'corr', 'znoise'});
    KeepTime.BadChannelDetection = toc(D);

    %%% Bad channel figure
    badchanFigDir = fullfile(opts.figpath, 'badchans', subDir);
    if ~exist(badchanFigDir, 'dir'), mkdir(badchanFigDir); end
    gedai.plotBadChannels(corr, znoise, EEG.chanlocs, ...
        fullfile(badchanFigDir, [fileID '_BadChannelTopoplot.png']));

    %%% Leadfield covariance matrix
    lfCOV = gedai.loadrefcov(opts.leadfielddir, p, EEG.nbchan, removed_channels);

    %%% Remove bad channels
    EEG = pop_select(EEG, 'nochannel', find(removed_channels));

    %%% Epochs to plot
    epochsToPlot = gedai.resolveEpochsToPlot(opts.epochstoplot, scoringDigits);

    %%% GEDAI runs
    for iRun = 1:numel(opts.runs)
        r        = opts.runs{iRun};
        savename = gedai.buildSaveName(r, EEG.srate);

        switch opts.runmode
            case 'StageSpecific'
                savename        = [opts.prefix 'StageSpecific_' savename];
                stageLogic      = {[-2], [-3], [0], [1]};
                refCOV_perStage = {lfCOV, lfCOV, lfCOV, lfCOV};
            case 'StateWise'
                savename        = [opts.prefix 'StateWise_' savename];
                stageLogic      = {[-2, -3], [0, 1]};
                refCOV_perStage = {lfCOV, lfCOV};
            case 'WholeNight'
                savename        = [opts.prefix 'WholeNight_' savename];
                stageLogic      = {[-3:1]};
                refCOV_perStage = {lfCOV};
        end

        %%% Drop stage groups absent from scoring
        presentStages   = unique(scoringDigits_NoN1);
        keep            = cellfun(@(s) any(ismember(s, presentStages)), stageLogic);
        stageLogic      = stageLogic(keep);
        refCOV_perStage = refCOV_perStage(keep);

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
                'GEDAIMode',                  r.GEDAIMode, ...
                'GEDAIEpochSize',             r.GEDAIEpochSize, ...
                'GEDAILowCutOffFreq',         r.GEDAILowCutOffFreq, ...
                'BBEpochSize',                r.GEDAIBroadbandEpochSize, ...
                'BroadbandOnly',              r.broadbandOnly, ...
                'GEDAIEnovaChannelThreshold', r.GEDAIEnovaChannelThreshold, ...
                'PercentileThreshold',        r.percentileThreshold, ...
                'BBMinThreshold',             r.BBMinThreshold, ...
                'ComputeSENSAI',              r.computeSENSAI, ...
                'ICAtype',                    r.ICAtype, ...
                'RefCOV',                     refCOV_perStage), ...
            gedaiDatFile, opts.refresh, {'EEGgedai', '', '', ''});
        KeepTime.GEDAI = toc(D);
        fprintf('GEDAI took %.2f min\n', KeepTime.GEDAI / 60)

        %%% JSON sidecar: GEDAI parameters + timing (only written on fresh runs)
        if isNewRun
            [~, bvBase] = fileparts(gedaiDatFile);
            sidecarjson(KeepTime, ...
                fullfile(gedaiRunDir, [bvBase '.json']), ...
                struct('GEDAIParameters', r));
        end

        %%% Evaluation figures
        if ~exist(gedaiFigDir, 'dir'), mkdir(gedaiFigDir); end
        run.eval_clean(EEG, EEGgedai, scoringDigits_NoN1, ...
            'EpochLength', opts.epochlength, 'WelchWindow', 4, ...
            'EpochsToPlot', epochsToPlot, 'refresh', opts.refresh, ...
            'SavePath', fullfile(gedaiFigDir, fileID));
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

    catch ME
        fprintf('[ERROR] %s: %s\n', fileID, ME.message);
        failures{end+1} = struct('fileID', fileID, 'message', ME.message, 'report', ME.getReport()); %#ok<AGROW>
    end
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
