function GEDAI_BIDS(BIDS, opts)
%GEDAI_BIDS  Run GEDAI artefact-rejection on a BIDS sleep EEG dataset.
%
%   GEDAI_BIDS(BIDS)
%   GEDAI_BIDS(BIDS, Name, Value, ...)
%
%   Required
%   --------
%   BIDS              bids.layout object.
%
%   Non-standard BIDS paths  (all optional)
%   ----------------------------------------
%   scoringpath       Directory containing sleep-scoring files (.json or .csv).
%                     Default: <BIDS root>/derivatives/scoring/scores/Manual_Checked
%   sfppath           BIDS root path, passed to the study-specific SFP
%                     resolver (dispatched via fileID on "PM" or "DROP").
%                     Default: <BIDS root>
%   leadfielddir      Root directory for brainstorm leadfields, organised as
%                     <sub>/<ses>/headmodel_surf_openmeeg.mat.
%                     Default: <BIDS root>/../Data_Analysis/Brainstorm_db/Leadfield_PM/data
%
%   Saving
%   ------
%   savepath          Output directory. Default: <BIDS root>/derivatives/GEDAI
%   refresh           Force re-run even if a cache file exists. Default: false.
%
%   EEG
%   ---
%   tasklabel         BIDS task label(s) to query. Default: {'Sleep','sleep'}.
%   recordinglabel    BIDS recording label to query. Default: '125Hz'.
%   noteegchannels    Channel indices to drop before processing.
%                     Default: [257:264].
%   targetsrate       Resample EEG to this rate (Hz). 0 = no resampling.
%                     Default: 125.
%   net               EEG net identifier passed to chans1020. Default: 'EGI256'.
%
%   Subject filter
%   --------------
%   subjectfilter     Cell array of subject ID strings; process only files
%                     whose fileID contains at least one entry.
%                     Default: {'hpmam006'}.
%
%   GEDAI
%   -----
%   runmode           'StageSpecific' (separate GEDAI per stage) or
%                     'WholeNight' (single run over all stages).
%                     Default: 'StageSpecific'.
%   epochlength       Sleep-epoch length in seconds. Default: 30.
%   runs              Cell array of GEDAI run-config structs. Default: one
%                     run with ICAtype='none' and canonical settings.
%   epochstoplot      Epoch indices included in diagnostic figures.
%                     Default: auto-selected (first + middle of each stage).
%
%   See dependancies.m.

arguments
    BIDS   % bids.layout object

    %--- Non-standard BIDS paths ---
    opts.scoringpath      char = fullfile(BIDS.pth, 'derivatives\scoring\scores\Manual_Checked')
    opts.sfppath          char = BIDS.pth
    opts.leadfielddir     char = fullfile(BIDS.pth, '..', 'Data_Analysis\Brainstorm_db\Leadfield_PM\data')

    %--- Saving ---
    opts.savepath         char = fullfile(BIDS.pth, 'derivatives', 'GEDAI')
    opts.refresh (1,1) logical = false

    %--- EEG ---
    opts.tasklabel               = {'Sleep', 'sleep'}
    opts.recordinglabel     char = '125Hz'
    opts.noteegchannels          = [257:264]
    opts.targetsrate (1,1) double = 125
    opts.net                char = 'EGI256'
    opts.removeDC  (1,1) logical = true;
    opts.removeLN  (1,1) logical = true;

    %--- Subject filter ---
    opts.subjectfilter      cell = {}

    %--- GEDAI ---
    opts.runmode {mustBeMember(opts.runmode, {'WholeNight', 'StageSpecific', 'StateWise'})} = 'StageSpecific'
    opts.epochlength (1,1) double = 30
    opts.runs                    = []
    opts.epochstoplot            = []
    opts.prefix             char = '';
end


%%% Initiate variables
KeepTime = [];
if isempty(opts.runs), opts.runs = gedai.defaultRuns(); end

%%% Query EEG files
filesEEG = bids.query(BIDS, 'data', 'extension', '.vhdr', ...
    'task', opts.tasklabel, 'recording', opts.recordinglabel);
if isempty(filesEEG); error('GEDAI_BIDS:noFiles', 'No matching EEG files found in BIDS layout.'); end

%%% Scoring files
scoringfiles = gedai.collectScoringFiles(opts.scoringpath);

%%% Loop over EEG files
for ifile = 1:numel(filesEEG)
    eegFile  = filesEEG{ifile};
    p        = bids.internal.parse_filename(eegFile);
    fileID   = strjoin(cellfun(@(k) [k '-' p.entities.(k)], fieldnames(p.entities), 'uni', 0), '_');

    %%% Subject filter
    if ~isempty(opts.subjectfilter) && ~contains(fileID, opts.subjectfilter)
        continue
    end
    fprintf('\n=== %s ===\n', fileID)
    fprintf('Save path → %s\n', opts.savepath)

    %%% Find matching scoring file (match on sub+ses only; recording label may differ)
    scoringFile = gedai.matchScoringFile(p.entities, scoringfiles);
    if isempty(scoringFile); error('GEDAI_BIDS:noScoring', 'No scoring file matched for %s.', fileID); end

    %%% Load sleep scoring
    fprintf('\nReading %s ...\n', scoringFile)
    scoringDigits = scoreloader(scoringFile);

    %%% Import EEG
    EEG = eeg_import(eegFile);

    %%% Correct scoring length if needed
    nEpochs = floor(EEG.pnts / (30 * EEG.srate));
    while numel(scoringDigits) > nEpochs; scoringDigits(end) = []; end

    %%% Drop non-EEG channels
    EEG = pop_select(EEG, 'nochannel', intersect(1:EEG.nbchan, opts.noteegchannels));

    %%% Read SFP file (dome-solved channel locations)
    if ~isempty(opts.sfppath)
        sfpFile = gedai.matchSfpFile(opts.sfppath, p.entities.sub, p.entities.ses);
        fprintf('\nReading %s ...\n', sfpFile)
        chanlocs = readlocs(sfpFile);
        chanlocs_reg = register_fiducials(chanlocs);
        EEG.chanlocs = chanlocs_reg(1:EEG.nbchan);
%         EEG.chanlocs = chanlocs(4 : EEG.nbchan + 3);
%     else
%         EEG.chanlocs = readlocs('C:\Postdoc\Code\exploratory-prep\locfiles\electrodes.tsv')
    end
EEG.chanlocs
    %%% Optional downsampling
    if opts.targetsrate > 0 && EEG.srate ~= opts.targetsrate
        fprintf('Resampling %d → %d Hz ...\n', EEG.srate, opts.targetsrate)
        EEG = pop_resample(EEG, opts.targetsrate);
    end
EEG.chanlocs
    %%% Build filters
    fprintf('Building filters (srate = %d Hz) ...\n', EEG.srate)
    EEG_DCFilter_NumDen = filterbank(EEG.srate, 'DC_RCSquareFilt');
    EEG_NotchFilt_IIR   = filterbank(EEG.srate, 'EEG_NotchFilt_IIR2');

    %%% DC removal
    if opts.removeDC
        D = tic; fprintf('\nDC removal ...\n')
        EEG.data = filtfilt(EEG_DCFilter_NumDen(1,:), EEG_DCFilter_NumDen(2,:), double(EEG.data'))';
        KeepTime.DCRemoval = toc(D);
    end
EEG.chanlocs
    %%% Notch filter
    if opts.removeLN
        D = tic; EEG.data = double(EEG.data);
        for ifilt = 1:numel(EEG_NotchFilt_IIR)
            fprintf('%d Hz notch ...\n', 50 * ifilt)
            EEG.data = filtfilt(EEG_NotchFilt_IIR{ifilt}, EEG.data')';
        end
        KeepTime.NotchFilter = toc(D);
    end
EEG.chanlocs
    %%% Average re-reference
    EEG.data = EEG.data - sum(EEG.data, 1) / (size(EEG.data, 1) + 1);
EEG.chanlocs
    %%% Replace isolated N1 epochs at stage boundaries with neighbour stage
    scoringDigits_NoN1 = gedai.killN1(scoringDigits);

    %%% Assign 10-20 labels
    EEG = chans1020(EEG, false, 'net', opts.net);

    %%% Bad channel detection
    D = tic;
    [removed_channels, corr, znoise] = smartcache( ...
        @() clean_channels(EEG, 0.7, 4, [], 0.5, 25), ...
        fullfile(opts.savepath, 'BadChannels', ['BadChans_' fileID '.mat']), ...
        false, {'', 'removed_channels', 'corr', 'znoise'});
    KeepTime.BadChannelDetection = toc(D);

    %%% Bad channel plot
    gedai.plotBadChannels(corr, znoise, EEG.chanlocs, fullfile(opts.savepath, 'BadChannels', [fileID '_BadChannelTopoplot.png']));

    %%% Leadfield covariance matrix
    lfCOV = gedai.loadrefcov(opts.leadfielddir, p, EEG.nbchan, removed_channels);

    %%% Remove bad channels
    EEG.urchanlocs = EEG.chanlocs;
    EEG = pop_select(EEG, 'nochannel', find(removed_channels));

    %%% Epochs to plot
    epochsToPlot = gedai.resolveEpochsToPlot(opts.epochstoplot, scoringDigits);

    %%% GEDAI runs
    for iRun = 1:numel(opts.runs)
        r        = opts.runs{iRun}; % GEDAI parameters
        savename = gedai.buildSaveName(r, EEG.srate);
        fprintf('Run %d/%d: %s\n', iRun, numel(opts.runs), savename)

        % Define GEDAI logic
        clear EEGgedai
        switch opts.runmode
            case 'StageSpecific'
                savename = [opts.prefix 'StageSpecific_' savename];
                stageLogic = {[-2], [-3], [0], [1]}; 
                refCOV_perStage = {lfCOV, lfCOV, lfCOV, lfCOV};

            case 'StateWise'
                savename = [opts.prefix 'StateWise_' savename];                
                stageLogic = {[-2, -3], [0, 1]}; 
                refCOV_perStage = {lfCOV, lfCOV};

            case 'WholeNight'
                savename = [opts.prefix 'StateWise_' savename];
                stageLogic = {[-3:1]}; 
                refCOV_perStage = {lfCOV};
        end
             
        % Run GEDAI
        [EEGgedai, ndxepochs, KeepTime] = smartcache( ...
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
            fullfile(opts.savepath, 'EEG', savename, [fileID '.mat']), ...
            opts.refresh, {'EEGgedai', '', 'ndxepochs', 'KeepTime'});
        fprintf('GEDAI took %.2f min\n', KeepTime.GEDAI / 60)

        % temporary stuff
%         EEG             = pop_interp(EEG, EEG.urchanlocs, 'spherical');
%         EEGgedai        = pop_interp(EEGgedai, EEG.urchanlocs, 'spherical');
%         chanlocs        = readlocs('C:\Postdoc\Code\exploratory-prep\locfiles\electrodes.tsv');  
%         chanlocs        = fixchanlocs(chanlocs);
%         EEG.chanlocs        = chanlocs(1:256);
%         EEGgedai.chanlocs   = chanlocs(1:256);
        EEGgedai.chanlocs = EEG.chanlocs;
        EEGgedai.urchanlocs = EEG.urchanlocs;
%         EEG.urchanlocs      = chanlocs(1:256);
%         EEGgedai.urchanlocs= chanlocs(1:256);
        savename = ['RegFiducials_' savename];
        epochsToPlot = [];
        EEG =eeg_checkset(EEG); EEGgedai=eeg_checkset(EEGgedai);
        run.eval_clean(EEG, EEGgedai, scoringDigits_NoN1(ndxepochs), ...
            'EpochLength', opts.epochlength, 'WelchWindow', 4, ...
            'EpochsToPlot', epochsToPlot, 'refresh', false, ...
            'SavePath', fullfile(opts.savepath, 'Figures', [savename], fileID, fileID))
        close all;


        %%% Evaluation plot
        run.eval_clean(EEG, EEGgedai, scoringDigits_NoN1(ndxepochs), ...
            'EpochLength', opts.epochlength, 'WelchWindow', 4, ...
            'EpochsToPlot', epochsToPlot, 'refresh', opts.refresh, ...
            'SavePath', fullfile(opts.savepath, 'Figures', [savename], fileID, fileID))
        close all;

        % Will be stand-alone soon
        if isfield(EEGgedai.etc, 'ic_classification')
            EEGgedai = ica.selectcomps(EEGgedai, 'ArtefactThreshold', 0.5, 'ManualQC', false);
            EEGgedai = pop_subcomp(EEGgedai, find(EEGgedai.reject.gcompreject), 0);

            run.eval_clean(EEG, EEGgedai, scoringDigits_NoN1(ndxepochs), ...
                'EpochLength', opts.epochlength, 'WelchWindow', 4, ...
                'EpochsToPlot', epochsToPlot, 'refresh', opts.refresh, ...
                'SavePath', fullfile(opts.savepath, 'Figures', savename, fileID, fileID))
            close all;
        end
    end
end

end % GEDAI_BIDS

