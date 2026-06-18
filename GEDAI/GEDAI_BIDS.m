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
%   sfppath           Root of the Data_collection tree, used to locate the
%                     dome-solved .sfp file for each subject/session.
%                     Default: <BIDS root>/../Data_collection
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
    opts.sfppath          char = fullfile(BIDS.pth, '..', 'Data_collection')
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

    %--- Subject filter ---
    opts.subjectfilter      cell = {}

    %--- GEDAI ---
    opts.runmode {mustBeMember(opts.runmode, {'WholeNight', 'StageSpecific'})} = 'StageSpecific'
    opts.epochlength (1,1) double = 30
    opts.runs                    = []
    opts.epochstoplot            = []
    opts.prefix             char = '';
end


%%% Initiate variables
KeepTime = [];
if isempty(opts.runs), opts.runs = defaultRuns(); end

%%% Query EEG files
filesEEG = bids.query(BIDS, 'data', 'extension', '.vhdr', ...
    'task', opts.tasklabel, 'recording', opts.recordinglabel);
if isempty(filesEEG); error('GEDAI_BIDS:noFiles', 'No matching EEG files found in BIDS layout.'); end

%%% Scoring files
scoringfiles = collectScoringFiles(opts.scoringpath);

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
    scoringFile = matchScoringFile(p.entities, scoringfiles);
    if isempty(scoringFile); error('GEDAI_BIDS:noScoring', 'No scoring file matched for %s.', fileID); end

    %%% Load sleep scoring
    fprintf('\nReading %s ...\n', scoringFile)
    scoringDigits = scoreloader(scoringFile);

    %%% Import EEG
    EEG = eeg_import(eegFile);

    %%% Drop non-EEG channels
    EEG = pop_select(EEG, 'nochannel', intersect(1:EEG.nbchan, opts.noteegchannels));

    %%% Read SFP file (dome-solved channel locations)
    if ~isempty(opts.sfppath)
        sfpFile = matchSfpFile(opts.sfppath, p.entities.sub, p.entities.ses);
        fprintf('\nReading %s ...\n', sfpFile)
        chanlocs = readlocs(sfpFile);
        EEG.chanlocs = chanlocs(4 : EEG.nbchan + 3);
    end

    %%% Optional downsampling
    if opts.targetsrate > 0 && EEG.srate ~= opts.targetsrate
        fprintf('Resampling %d → %d Hz ...\n', EEG.srate, opts.targetsrate)
        EEG = pop_resample(EEG, opts.targetsrate);
    end

    %%% Build filters
    fprintf('Building filters (srate = %d Hz) ...\n', EEG.srate)
    EEG_DCFilter_NumDen = filterbank(EEG.srate, 'DC_RCSquareFilt');
    EEG_NotchFilt_IIR   = filterbank(EEG.srate, 'EEG_NotchFilt_IIR2');

    %%% DC removal
    D = tic; fprintf('\nDC removal ...\n')
    EEG.data = filtfilt(EEG_DCFilter_NumDen(1,:), EEG_DCFilter_NumDen(2,:), double(EEG.data'))';
    KeepTime.DCRemoval = toc(D);

    %%% Notch filter
    D = tic; EEG.data = double(EEG.data);
    for ifilt = 1:numel(EEG_NotchFilt_IIR)
        fprintf('%d Hz notch ...\n', 50 * ifilt)
        EEG.data = filtfilt(EEG_NotchFilt_IIR{ifilt}, EEG.data')';
    end
    KeepTime.NotchFilter = toc(D);

    %%% Average re-reference
    EEG.data = EEG.data - sum(EEG.data, 1) / (size(EEG.data, 1) + 1);

    %%% Replace isolated N1 epochs at stage boundaries with neighbour stage
    scoringDigits_NoN1 = killN1(scoringDigits);

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
    plotBadChannels(corr, znoise, EEG.chanlocs, fullfile(opts.savepath, 'BadChannels', 'BadChannelTopoplot.png'));

    %%% Leadfield covariance matrix
    lfCOV = loadrefcov(opts.leadfielddir, p, EEG.nbchan, removed_channels);

    %%% Remove bad channels
    EEG = pop_select(EEG, 'nochannel', find(removed_channels));

    %%% Epochs to plot
    epochsToPlot = resolveEpochsToPlot(opts.epochstoplot, scoringDigits);

    %%% GEDAI runs
    for iRun = 1:numel(opts.runs)
        r        = opts.runs{iRun};
        savename = buildSaveName(r, EEG.srate);
        fprintf('Run %d/%d: %s\n', iRun, numel(opts.runs), savename)

        clear EEGgedai
        switch opts.runmode
            case 'StageSpecific'
                savename = [opts.prefix 'StageSpecific_' savename];
                [EEGgedai, ndxepochs, KeepTime] = smartcache( ...
                    @() run.GEDAI_StageSpecific(EEG, scoringDigits_NoN1, ...
                        {[-2], [-3], [0], [1]}, KeepTime, ...
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
                        'RefCOV',                     {lfCOV, lfCOV, lfCOV, lfCOV}), ...
                    fullfile(opts.savepath, 'EEG', [savename '_' fileID '.mat']), ...
                    opts.refresh, {'EEGgedai', '', 'ndxepochs', 'KeepTime'});

            case 'WholeNight'
                savename = [opts.prefix 'WholeNight_' savename];                
                [EEGgedai, ndxepochs, KeepTime] = smartcache( ...
                    @() run.GEDAI_StageSpecific(EEG, scoringDigits_NoN1, ...
                        {[-3:1]}, KeepTime, ...
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
                        'RefCOV',                     {lfCOV}), ...
                    fullfile(opts.savepath, 'EEG', [savename '_' fileID '.mat']), ...
                    opts.refresh, {'EEGgedai', '', 'ndxepochs', 'KeepTime'});
        end
        fprintf('GEDAI took %.2f min\n', KeepTime.GEDAI / 60)

        %%% Evaluation plot
        run.eval_clean(EEG, EEGgedai, scoringDigits_NoN1(ndxepochs), ...
            'EpochLength', opts.epochlength, 'WelchWindow', 4, ...
            'EpochsToPlot', epochsToPlot, 'refresh', opts.refresh, ...
            'SavePath', fullfile(opts.savepath, 'Figures', [savename], fileID, fileID))
        close all;

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


%% =========================================================================
%  LOCAL HELPERS
% ==========================================================================

function runs = defaultRuns()
    runs = {struct( ...
        'GEDAIMode',                  'auto', ...
        'GEDAILowCutOffFreq',         0.1, ...
        'GEDAIEpochSize',             12, ...
        'GEDAIBroadbandEpochSize',    10, ...
        'boost1',                     1, ...
        'boost2',                     1, ...
        'broadbandOnly',              false, ...
        'percentileThreshold',        98, ...
        'WeightKC',                   0, ...
        'BBMinThreshold',             -2, ...
        'computeSENSAI',              false, ...
        'GEDAIEnovaChannelThreshold', Inf, ...
        'ICAtype',                    'none')};
end


function sfp = matchSfpFile(CollectionRoot, SubjectName, SessName)
% Map a BIDS (subject, session) to its .sfp in the Data_collection tree.
%
% [VERIFY] Inferred from two (non-matching) examples:
%     sub-hpmam003  ->  H003_PM_AM           (trailing number fills the skeleton)
%     ses-S1        ->  S1
%     .../H###_PM_AM/SA_stim/S#/GPS/H###_PM_AM_GPS_S#_coordinates.sfp
% Edit this single function if your naming differs.
    lab    = erase(SubjectName, 'sub-');
    num    = regexp(lab, '\d+$', 'match', 'once');
    if isempty(num), sfp = ''; return; end
    subColl = sprintf('H%s_PM_AM', num);
    ses     = erase(SessName, 'ses-');
    gpsDir  = fullfile(CollectionRoot, subColl, 'SA_stim', ses, 'GPS');
    sfp     = fullfile(gpsDir, sprintf('%s_GPS_%s_coordinates.sfp', subColl, ses));
    if ~isfile(sfp)
        d = dir(fullfile(gpsDir, '*coordinates*.sfp'));
        if isempty(d), d = dir(fullfile(gpsDir, '*.sfp')); end
        if ~isempty(d), sfp = fullfile(d(1).folder, d(1).name); end
    end
end


function scoringFile = matchScoringFile(entities, allScoring)
% Match by sub + ses only (recording label may differ between EEG and scoring).
    scoringFile = '';
    if isempty(allScoring), return; end
    subs = regexp(allScoring, '(?<=sub-)[^_]+', 'match', 'once');
    sess = regexp(allScoring, '(?<=ses-)[^_]+', 'match', 'once');
    idx  = find(strcmp(subs, entities.sub) & strcmp(sess, entities.ses), 1);
    if ~isempty(idx), scoringFile = allScoring{idx}; end
end


function digits = killN1(scoringDigits)
% Replace isolated N1 epochs at stage boundaries with their neighbour's stage.
    digits = scoringDigits;
    while any(digits == -1)
        firstN1 = find([0,  diff(digits == -1)] ==  1);
        lastN1  = find([diff(digits == -1),  0] == -1);
        digits(firstN1) = digits(firstN1 - 1);
        digits(lastN1)  = digits(lastN1  + 1);
    end
end


function epochsToPlot = resolveEpochsToPlot(requested, scoringDigits)
    epochsToPlot = requested(:)';
    for score = -3:1
        idx = find(scoringDigits == score);
        if isempty(idx), continue; end
        epochsToPlot = [epochsToPlot, idx(1), idx(round(end/2)), idx(round(end/3)), idx(round(end/4)), idx(round(end/5))]; %#ok<AGROW>
        if score == -2 && numel(idx) >= 20
            epochsToPlot = [epochsToPlot, idx(5:5:20)]; %#ok<AGROW>
        end
    end
    epochsToPlot = unique(epochsToPlot);
end


function lfCOV = loadrefcov(leadfielddir, p, nbchan, removed_channels)

    if isfolder(leadfielddir)
        lfFile = fullfile(leadfielddir, ['sub-' p.entities.sub], ['ses-' p.entities.ses], 'headmodel_surf_openmeeg.mat');
    end
    bstorm    = load(lfFile);
    goodChans = setdiff(1:nbchan, find(removed_channels));
    B         = bstorm.Gain(goodChans, :);
    B         = B - sum(B, 1) / (size(B, 1) + 1);
    lfCOV     = B * B';
end


function name = buildSaveName(r, srate)
    name = sprintf( ...
        'GEDAIBIDS_ICA%s_BBonly%d_LowCutOff%d_BBEpochSize%d_%s_PrcThresh%d_BBMinThresh%d_EnovaChanThresh%d_SENSAI%d_KC%d_b1x%d_b2x%d_Srate%d', ...
        r.ICAtype, r.broadbandOnly, r.GEDAILowCutOffFreq * 10, ...
        r.GEDAIBroadbandEpochSize, r.GEDAIMode, r.percentileThreshold, ...
        r.BBMinThreshold, r.GEDAIEnovaChannelThreshold, r.computeSENSAI, ...
        r.WeightKC * 100, round(r.boost1 * 10), round(r.boost2 * 10), srate);
end


function files = collectScoringFiles(scoringBase)
% JSON preferred, CSV fallback.
    f = dir(fullfile(scoringBase, '**', '*.json'));
    if isempty(f), f = dir(fullfile(scoringBase, '**', '*.csv')); end
    if isempty(f), files = {}; return; end
    files = fullfile({f.folder}, {f.name})';
end


function plotBadChannels(corr, znoise, chanlocs, savefile)
    lowcorrprop = sum(corr < .7, 2) ./ size(corr, 2);
    figure();
    tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    nexttile();
    topoplot(lowcorrprop, chanlocs, 'numcontour', 0, 'emarker2', {find(lowcorrprop > .5), '.', 'r', 10});
    colorbar(); caxis([0 .8]); title({'Prop. of recording', 'with low correlation'});
    nexttile();
    topoplot(znoise, chanlocs, 'numcontour', 0, 'emarker2', {find(znoise > 4), '.', 'r', 10});
    colorbar(); caxis([0 4]); title('Line noise');
    colormap('gray');
    set(gcf, 'Color', 'w', 'Units', 'centimeters', 'Position', [2 2 20 10]);
    print(gcf, savefile, '-dpng', '-r150');
    close
end
