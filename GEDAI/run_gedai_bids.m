function run_gedai_bids(BIDS, opts)
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
%                     Default: <BIDS root>/derivatives/preprocessing/A_filtered
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
%                     Default: <BIDS root>/derivatives/preprocessing/GEDAI
%   figpath           Root for all figures.
%                     Default: <BIDS root>/derivatives/preprocessing/figures
%   geddesc           desc label for GEDAI BrainVision output files.
%                     Default: 'filtGEDAI'
%   refresh           Force re-run even if a cache file exists. Default: false.
%
%   EEG
%   ---
%   tasklabel         BIDS task label(s) to query. Default: {'Sleep','sleep'}.
%   recordinglabel    BIDS recording label to query. Default: '125Hz'.
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

    %--- Input paths ---
    opts.inputpath     char = fullfile(BIDS.pth, 'derivatives', 'preprocessing')
    opts.inputdesc         char = 'filt'
    opts.scoringpath      char = fullfile(BIDS.pth, 'derivatives', 'scoring', 'scores', 'Manual_Checked')
    opts.sfppath          char = BIDS.pth
    opts.leadfielddir     char = fullfile(BIDS.pth, '..', 'Data_Analysis', 'Brainstorm_db', 'Leadfield_PM', 'data')

    %--- Output paths ---
    opts.savepath         char = fullfile(BIDS.pth, 'derivatives', 'preprocessing')
    opts.figpath          char = fullfile(BIDS.pth, 'derivatives', 'preprocessing', 'figures')
    opts.geddesc          char = 'filtGEDAI'
    opts.refresh (1,1) logical = false

    %--- EEG ---
    opts.tasklabel                      = {'Sleep', 'sleep'}
    opts.recordinglabel    char         = '125Hz'
    opts.noteegchannels    (1,:) double = 257:264
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

KeepTime = struct();
if isempty(opts.runs), opts.runs = gedai.defaultRuns(); end

%%% Query EEG files from BIDS (used for entity extraction and subject iteration)
filesEEG = bids.query(BIDS, 'data', 'extension', '.vhdr', ...
    'task', opts.tasklabel, 'acq', opts.recordinglabel);
if isempty(filesEEG)
    error('run_gedai_bids:noFiles', 'No matching EEG files found in BIDS layout.');
end

%%% Scoring files
scoringfiles = gedai.collectScoringFiles(opts.scoringpath);

%%% Loop over EEG files
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

    %%% Resolve filtered input file
    filtFile = fullfile(opts.inputpath, subDir, [fileID '_desc-' opts.inputdesc '_eeg.vhdr']);
    if ~isfile(filtFile)
        fprintf('[skip] %s: filtered file not found (%s)\n', fileID, filtFile)
        continue
    end
    fprintf('Input  → %s\n', filtFile)
    fprintf('Output → %s\n', opts.savepath)

    %%% Find matching scoring file
    scoringFile = gedai.matchScoringFile(p.entities, scoringfiles);
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
    else
        continue
    end

    %%% Average re-reference
    EEG.data = EEG.data - sum(EEG.data, 1) / (size(EEG.data, 1) + 1);

    %%% Replace isolated N1 epochs at stage boundaries
    scoringDigits_NoN1 = gedai.killN1(scoringDigits);

    %%% Assign 10-20 labels
    EEG = chans1020(EEG, false, 'net', opts.net);

    %%% Bad channel detection
    badchanDir    = fullfile(opts.savepath, subDir, 'BadChannels');
    badchanFigDir = fullfile(opts.figpath, 'badchans', subDir);

    D = tic;
    [removed_channels, corr, znoise] = smartcache( ...
        @() clean_channels(EEG, 0.7, 4, [], 0.5, 25), ...
        fullfile(badchanDir, ['BadChans_' fileID '.mat']), ...
        false, {'', 'removed_channels', 'corr', 'znoise'});
    KeepTime.BadChannelDetection = toc(D);

    %%% Bad channel figure
    if ~exist(badchanFigDir, 'dir'), mkdir(badchanFigDir); end
    gedai.plotBadChannels(corr, znoise, EEG.chanlocs, ...
        fullfile(badchanFigDir, [fileID '_BadChannelTopoplot.png']));

    %%% Leadfield covariance matrix
    lfCOV = gedai.loadrefcov(opts.leadfielddir, p, EEG.nbchan, removed_channels);

    %%% Remove bad channels
    EEG.urchanlocs = EEG.chanlocs;
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
        fprintf('Run %d/%d: %s\n', iRun, numel(opts.runs), savename)

        %%% Paths for this run
        gedaiRunDir  = fullfile(opts.savepath, subDir, savename);
        gedaiFigDir  = fullfile(opts.figpath, 'gedai', savename, fileID);
        gedaiDatFile = fullfile(gedaiRunDir, [fileID '_desc-' opts.geddesc '_eeg.dat']);

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
end

end % run_gedai_bids
