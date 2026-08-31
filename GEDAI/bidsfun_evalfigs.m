function failures = bidsfun_evalfigs(BIDS, opts)
% BIDSFUN_EVALFIGS  Evaluate processing quality by comparing a "before" and "after" EEG file.
%
%   bidsfun_evalfigs(BIDS)
%   bidsfun_evalfigs(BIDS, Name, Value, ...)
%
%   Runs run.eval_clean on a paired before/after EEG file to visualise the effect of
%   whichever processing stage produced the "after" file.
%
%   Required
%   --------
%   BIDS              bids.layout object.
%
%   Input paths
%   -----------
%   derivpath         Root of the derivative files being compared.
%                     Default: <BIDS root>/derivatives/prep-ged
%   afterdesc         desc label of the "after" file, and of the figure folder.
%                     Default: 'filt'
%   beforedesc        desc label of the "before" file. '' (default) uses the raw BIDS
%                     recording; set it to compare two derivatives - e.g. beforedesc the
%                     filtered data and afterdesc the GEDAI output, which reproduces what
%                     bidsfun_gedai plots internally.
%   avgrefbefore      average-reference the "before" recording before plotting. Match the
%                     stage being reproduced: bidsfun_gedai re-references before GEDAI,
%                     bidsfun_hp_zap_cleanline does not. Default: false
%   avgrefafter       average-reference the "after" recording before plotting. Default: false
%   scoringpath       Directory containing sleep-scoring files.
%                     Default: <BIDS root>/derivatives/scoring/scores/Manual_Checked
%   sfppath           Path passed to the SFP resolver. Default: <BIDS root>
%
%   Output paths
%   ------------
%   figpath           Root for evaluation figures.
%                     Default: <BIDS root>/derivatives/prep-ged/figures
%   refresh           Force re-run even if cached files exist. Default: false.
%
%   EEG
%   ---
%   tasklabel         BIDS task label(s) to query. Default: {'Sleep','sleep'}.
%   acqlabel    BIDS acq label to query. Default: '125Hz'.
%   noteegchannels    Channel indices to drop from raw EEG. Default: 257:300.
%   net               EEG net identifier passed to chans1020. Default: 'EGI256'.
%
%   Subject filter
%   --------------
%   subjectfilter     Cell array of subject ID strings; {} = all subjects.
%   sessionfilter     Cell array of session ID strings; {} = all sessions.
%
%   Epoch
%   -----
%   epochlength       Epoch duration in seconds. Default: 30.
%
%   Plots (forwarded to run.eval_clean)
%   ------------------------------------
%   All default true except PlotTopoBandPower, PlotEpochOverlay,
%   PlotExponentByStage, PlotSlopesTimecourse, which default false.
%   See run.eval_clean for what each plot shows.

arguments
    BIDS

    %--- Derivative folder ---
    opts.derivfolder      char = 'prep-ged'

    %--- Input paths ---
    opts.derivpath        char = ''
    opts.afterdesc        char = ''
    opts.beforedesc       char = ''
    opts.beforefileext    char = '.set'
    opts.afterfileext     char = '.set'
    opts.avgrefbefore (1,1) logical = false
    opts.avgrefafter  (1,1) logical = false
    opts.scoringpath      char = []
    opts.sfppath          char = BIDS.pth

    %--- Output paths ---
    opts.figpath          char = ''
    opts.refresh (1,1)    logical = false

    %--- EEG ---
    opts.tasklabel                      = {'Sleep', 'sleep'}
    opts.acqlabel   char         = '125Hz'
    opts.noteegchannels   (1,:) double = 257:300
    opts.net              char         = 'EGI256'
    opts.targetchannelcount (1,1) double = 256

    %--- Subject filter ---
    opts.subjectfilter    cell          = {}
    opts.sessionfilter    cell          = {}

    %--- Epoch ---
    opts.epochlength (1,1) double = 30
    opts.epochstoplot             = []
    

    %--- Plots (forwarded to run.eval_clean) ---
    opts.PlotCharacteristics  (1,1) logical = true
    opts.PlotPsdPerStage      (1,1) logical = true
    opts.PlotPsdPerStageChans (1,1) logical = true
    opts.PlotPsdOverview      (1,1) logical = true
    opts.PlotTopoBandPower    (1,1) logical = true
    opts.PlotTopoBandStage    (1,1) logical = true
    opts.PlotEpochOverlay     (1,1) logical = false
    opts.PlotTimefreq         (1,1) logical = true
    opts.PlotExponentByStage  (1,1) logical = false
    opts.PlotSlopesTimecourse (1,1) logical = false
end

fprintf('\n=== Running bidsfun_evalfigs ===\n');

if isempty(opts.derivpath), opts.derivpath = fullfile(BIDS.pth, 'derivatives', opts.derivfolder); end
if isempty(opts.figpath),      opts.figpath      = fullfile(BIDS.pth, 'derivatives', opts.derivfolder, 'figures'); end

%%% Query raw EEG files from BIDS
filesEEG = bids.query(BIDS, 'data', 'extension', '.vhdr', ...
    'task', opts.tasklabel, 'acq', opts.acqlabel);
if isempty(filesEEG)
    error('bidsfun_evalfigs:noFiles', 'No matching EEG files found in BIDS layout.');
end

%%% Scoring files
if ~isempty(opts.scoringpath)
    scoringfiles = gedai.collectScoringFiles(opts.scoringpath);
end

%%% Loop over EEG files
failures = {};
for ifile = 1:numel(filesEEG)
    rawFile = filesEEG{ifile};
    p       = bids.internal.parse_filename(rawFile);
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

    %%% Skip if already processed
    pngFile = fullfile(opts.figpath, ['desc-' opts.afterdesc], subDir, [fileID '_psd_per_stage.png']);
    if ~opts.refresh && isfile(pngFile)
        fprintf('[skip] output exists: %s\n', pngFile)
        continue
    end

    try

        %%% Resolve "after" file
        afterFile = fullfile(opts.derivpath, subDir, [fileID '_desc-' opts.afterdesc '_eeg' opts.afterfileext]);
        if ~isfile(afterFile)
            fprintf('[skip] after file not found: %s\n', afterFile)
            continue
        end
        fprintf('Raw    → %s\nAfter  → %s\n', rawFile, afterFile)

        %%% Find and load scoring
        if isempty(opts.scoringpath)
            scoringpath = fullfile(BIDS.pth, subDir);
            scoringfiles = gedai.collectScoringFiles(scoringpath);
        end
        scoringFile = gedai.matchScoringFile(p.entities, scoringfiles);
        if isempty(scoringFile)
            error('bidsfun_evalfigs:noScoring', 'No scoring file matched for %s.', fileID);
        end
        fprintf('Scoring → %s\n', scoringFile)
        scoringDigits = scoreloader(scoringFile);

        %%% Import raw EEG
        %%% The "before" recording: the raw BIDS file by default, or another derivative when
        %%% beforedesc is set - which is what lets this stage produce the GEDAI comparison
        %%% (before = the filtered input, after = the GEDAI output) and not only raw-vs-filtered.
        if isempty(opts.beforedesc)
            fprintf('Importing raw EEG ...\n')
            EEGraw = fast_eeg_import(rawFile);
        else
            beforeFile = fullfile(opts.derivpath, subDir, ...
                [fileID '_desc-' opts.beforedesc '_eeg' opts.beforefileext]);
            if ~isfile(beforeFile)
                fprintf('[skip] before file not found: %s\n', beforeFile)
                continue
            end
            fprintf('Before → %s\n', beforeFile)
            EEGraw = fast_eeg_import(beforeFile);
        end

        %%% Correct scoring length if needed
        nEpochs = floor(EEGraw.pnts / (opts.epochlength * EEGraw.srate));
        while numel(scoringDigits) > nEpochs; scoringDigits(end) = []; end


        %%% Channel locations. Coordinates come from the file when it carries them - a
        %%% derivative .set does, and after bad channel removal channel k is no longer the
        %%% k-th SFP entry, so gedai.assignChanlocs leaves those alone and only falls back
        %%% to the SFP/ercp sources when they are missing.
        if opts.targetchannelcount < EEGraw.nbchan
            EEGraw = pop_select(EEGraw, 'nochannel', intersect(1:EEGraw.nbchan, opts.noteegchannels));
        end
        EEGraw = gedai.assignChanlocs(EEGraw, BIDS, opts.sfppath, rawFile, p, fileID);
        
        %%% Import "after" EEG
        fprintf('Importing after EEG ...\n')
        EEGafter = fast_eeg_import(afterFile);
        if opts.targetchannelcount < EEGafter.nbchan        
            EEGafter = pop_select(EEGafter, 'nochannel', intersect(1:EEGafter.nbchan, opts.noteegchannels));
        end
        EEGafter = gedai.assignChanlocs(EEGafter, BIDS, opts.sfppath, rawFile, p, fileID);

        % %%% Removed channels
        % removed_channels = true(numel(EEGafter.urchanlocs), 1);
        % removed_channels([EEGafter.chanlocs.urchan]) = false;        
        % EEGraw = pop_select(EEGraw, 'nochannel', intersect(1:EEGraw.nbchan, find(removed_channels)));

        %%% Average reference, when asked, per recording. bidsfun_gedai re-references before
        %%% it runs GEDAI, so its figures are of average-referenced data; set avgrefafter (and
        %%% avgrefbefore, if the "before" file is itself a re-referenced derivative) to match
        %%% when reproducing that comparison, or the two stages' figures sit on different
        %%% baselines. Same formula as there - the +1 counts the implicit reference channel.
        if opts.avgrefbefore
            fprintf('Average-referencing before recording ...\n')
            EEGraw.data = EEGraw.data - sum(EEGraw.data, 1) / (size(EEGraw.data, 1) + 1);
            % EEGraw.data = EEGraw.data - sum(EEGraw.data(~removed_channels, :), 1) / (size(EEGraw.data, 1) + 1 - sum(removed_channels));
        end
        if opts.avgrefafter
            fprintf('Average-referencing after recording ...\n')
            EEGafter.data = EEGafter.data - sum(EEGafter.data, 1) / (size(EEGafter.data, 1) + 1);
        end

        %%% Epochs to plot
        epochsToPlot = gedai.resolveEpochsToPlot(opts.epochstoplot, scoringDigits);

        %%% Evaluate
        figDir = fullfile(opts.figpath, ['desc-' opts.afterdesc], subDir);
        if ~exist(figDir, 'dir'), mkdir(figDir); end

        run.eval_clean(EEGraw, EEGafter, scoringDigits, ...
            'EpochLength', opts.epochlength, 'WelchWindow', 4, ...
            'SavePath', fullfile(figDir, fileID), ...
            'refresh', opts.refresh, ...
            'net', opts.net, ...
            'EpochsToPlot', epochsToPlot, ...
            'PlotCharacteristics', opts.PlotCharacteristics, ...
            'PlotPsdPerStage', opts.PlotPsdPerStage, ...
            'PlotPsdPerStageChans', opts.PlotPsdPerStageChans, ...
            'PlotPsdOverview', opts.PlotPsdOverview, ...
            'PlotTopoBandPower', opts.PlotTopoBandPower, ...
            'PlotTopoBandStage', opts.PlotTopoBandStage, ...
            'PlotEpochOverlay', opts.PlotEpochOverlay, ...
            'PlotTimefreq', opts.PlotTimefreq, ...
            'PlotExponentByStage', opts.PlotExponentByStage, ...
            'PlotSlopesTimecourse', opts.PlotSlopesTimecourse);
        close all;

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
    if ~exist(opts.figpath, 'dir'), mkdir(opts.figpath); end
    fid = fopen(fullfile(opts.figpath, 'failed_files_evalfilt.json'), 'w');
    fprintf(fid, '%s', jsonencode([failures{:}], 'PrettyPrint', true));
    fclose(fid);
end
end
