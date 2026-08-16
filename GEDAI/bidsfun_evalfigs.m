function failures = bidsfun_evalfigs(BIDS, opts)
% BIDSFUN_EVALFIGS  Evaluate filtering quality by comparing raw and filtered EEG.
%
%   bidsfun_evalfigs(BIDS)
%   bidsfun_evalfigs(BIDS, Name, Value, ...)
%
%   Runs run.eval_clean on paired raw and desc-filt EEG files to
%   visualise the effect of the filtering pipeline.
%
%   Required
%   --------
%   BIDS              bids.layout object.
%
%   Input paths
%   -----------
%   filteredpath      Root of the filtered derivatives.
%                     Default: <BIDS root>/derivatives/prep-ged
%   filtdesc          desc label of the "after" file, and of the figure folder.
%                     Default: 'filt'
%   beforedesc        desc label of the "before" file. '' (default) uses the raw BIDS
%                     recording; set it to compare two derivatives - e.g. beforedesc the
%                     filtered data and filtdesc the GEDAI output, which reproduces what
%                     bidsfun_run_gedai plots internally.
%   avgref            average-reference both recordings before plotting. Match the stage
%                     being reproduced: bidsfun_run_gedai re-references before GEDAI,
%                     bidsfun_hp_zap_cleanline does not. Default: false
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
    opts.filteredpath     char = ''
    opts.filtdesc         char = 'filt'
    opts.beforedesc       char = ''
    opts.filtfileext      char = '.set'
    opts.avgref (1,1)     logical = false
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

    %--- Subject filter ---
    opts.subjectfilter    cell          = {}
    opts.sessionfilter    cell          = {}

    %--- Epoch ---
    opts.epochlength (1,1) double = 30

    %--- Plots (forwarded to run.eval_clean) ---
    opts.PlotCharacteristics  (1,1) logical = true
    opts.PlotPsdPerStage      (1,1) logical = true
    opts.PlotPsdPerStageChans (1,1) logical = true
    opts.PlotPsdOverview      (1,1) logical = true
    opts.PlotTopoBandPower    (1,1) logical = false
    opts.PlotTopoBandStage    (1,1) logical = true
    opts.PlotEpochOverlay     (1,1) logical = false
    opts.PlotTimefreq         (1,1) logical = true
    opts.PlotExponentByStage  (1,1) logical = false
    opts.PlotSlopesTimecourse (1,1) logical = false
end

if isempty(opts.filteredpath), opts.filteredpath = fullfile(BIDS.pth, 'derivatives', opts.derivfolder); end
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
    pngFile = fullfile(opts.figpath, ['desc-' opts.filtdesc], subDir, [fileID '_psd_per_stage.png']);
    if ~opts.refresh && isfile(pngFile)
        fprintf('[skip] output exists: %s\n', pngFile)
        continue
    end

    try

    %%% Resolve filtered file
    filtFile = fullfile(opts.filteredpath, subDir, [fileID '_desc-' opts.filtdesc '_eeg' opts.filtfileext]);
    if ~isfile(filtFile)
        fprintf('[skip] filtered file not found: %s\n', filtFile)
        continue
    end
    fprintf('Raw    → %s\nFilt   → %s\n', rawFile, filtFile)

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
        EEGraw = eeg_import(rawFile);
    else
        beforeFile = fullfile(opts.filteredpath, subDir, ...
            [fileID '_desc-' opts.beforedesc '_eeg' opts.filtfileext]);
        if ~isfile(beforeFile)
            fprintf('[skip] before file not found: %s\n', beforeFile)
            continue
        end
        fprintf('Before → %s\n', beforeFile)
        EEGraw = eeg_import(beforeFile);
    end

    %%% Correct scoring length if needed
    nEpochs = floor(EEGraw.pnts / (opts.epochlength * EEGraw.srate));
    while numel(scoringDigits) > nEpochs; scoringDigits(end) = []; end
    

    %%% Load SFP
    if strcmpi(BIDS.description.Name, {'ercp'})
        chanfile = fullfile(fileparts(rawFile), [fileID, '_channels.tsv']);
        elecfile = fullfile(fileparts(rawFile), ['sub-' p.entities.sub, '_ses-' p.entities.ses, '_electrodes.tsv']);
        [EEGraw, channelData, elecData] = bids_importchanlocs(EEGraw, chanfile, elecfile);
        chanlocs = EEGraw.chanlocs;
    else
        sfpFile      = gedai.matchSfpFile(opts.sfppath, p.entities.sub, p.entities.ses);
        fprintf('SFP     → %s\n', sfpFile)
        chanlocs = register_fiducials(readlocs(sfpFile));
    end
        
    %%% Extract channels. Coordinates come from the file when it carries them - a
    %%% derivative .set does, and after bad channel removal channel k is no longer the
    %%% k-th SFP entry, so overwriting would mislabel every channel. Otherwise fall back
    %%% to the SFP, which assumes the full montage in SFP order.
    EEGraw = pop_select(EEGraw, 'nochannel', intersect(1:EEGraw.nbchan, opts.noteegchannels));
    if ~isfield(EEGraw, 'chanlocs') || isempty(EEGraw.chanlocs) || ...
            ~isfield(EEGraw.chanlocs, 'X') || isempty([EEGraw.chanlocs.X])
        EEGraw.chanlocs   = chanlocs(1:EEGraw.nbchan);
        EEGraw.urchanlocs = EEGraw.chanlocs;
    end

    %%% Import filtered EEG
    fprintf('Importing filtered EEG ...\n')
    EEGfilt = eeg_import(filtFile);
    EEGfilt = pop_select(EEGfilt, 'nochannel', intersect(1:EEGfilt.nbchan, opts.noteegchannels));

    %%% Channel locations come from the .set itself. bidsfun_hp_zap_cleanline stored the
    %%% pre-removal montage in urchanlocs, and run.eval_clean interpolates the missing
    %%% channels back from it. Do not overwrite either here: after the removal channel k
    %%% is no longer the k-th SFP entry, so chanlocs(1:nbchan) would mislabel every
    %%% channel, and overwriting urchanlocs would destroy the record of what to restore.
    if ~isfield(EEGfilt, 'chanlocs') || isempty(EEGfilt.chanlocs) || ...
            ~isfield(EEGfilt.chanlocs, 'X') || isempty([EEGfilt.chanlocs.X])
        %%% Legacy input with no stored coordinates: assume the full montage, SFP order
        EEGfilt.chanlocs   = chanlocs(1:EEGfilt.nbchan);
        EEGfilt.urchanlocs = EEGfilt.chanlocs;
    end

    %%% Average reference, when asked. bidsfun_run_gedai re-references before it runs
    %%% GEDAI, so its figures are of average-referenced data; set avgref to match when
    %%% reproducing that comparison, or the two stages' figures sit on different
    %%% baselines. Same formula as there - the +1 counts the implicit reference channel.
    if opts.avgref
        fprintf('Average-referencing both recordings ...\n')
        EEGraw.data  = EEGraw.data  - sum(EEGraw.data,  1) / (size(EEGraw.data,  1) + 1);
        EEGfilt.data = EEGfilt.data - sum(EEGfilt.data, 1) / (size(EEGfilt.data, 1) + 1);
    end

    %%% Evaluate
    figDir = fullfile(opts.figpath, ['desc-' opts.filtdesc], subDir);
    if ~exist(figDir, 'dir'), mkdir(figDir); end

    run.eval_clean(EEGraw, EEGfilt, scoringDigits, ...        
        'EpochLength', opts.epochlength, 'WelchWindow', 4, ...
        'SavePath', fullfile(figDir, fileID), ...
        'refresh', opts.refresh, ...
        'net', opts.net, ...
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
