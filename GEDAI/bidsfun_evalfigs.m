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
%   Parallel machines
%   -----------------
%   Several machines can be pointed at the same BIDS root and the same figpath at once.
%   Each recording is claimed with a lock file before any work starts, so a second
%   machine walking the same list skips whatever the first is busy with instead of
%   duplicating it. The claim is released when the recording finishes, when it errors,
%   and on Ctrl-C, so a failure never parks a recording permanently. Claims are keyed by
%   afterdesc, so this stage never blocks a machine evaluating a different desc - or
%   running a different stage - on the same recording. See claimFile.
%
%   uselocks          Claim each recording before processing it. Default: true.
%                     Set false for a single-machine run over a local figpath.
%   lockpath          Directory holding the claim files. Default: <figpath>/.locks
%   lockstalemin      Minutes after which a claim whose heartbeat stopped is taken
%                     over - the escape hatch for a machine that crashed or was
%                     rebooted mid-recording. Keep it well above the longest
%                     plausible single-recording runtime. Default: 360 (6 h).
%   lockheartbeatmin  Minutes between heartbeat writes on a held claim. Default: 5.
%
%   Parallel pool
%   -------------
%   poolworkers       Workers for the FOOOF pool. 'auto' (default) sizes it from the
%                     machine's available memory and core count via gedai.autoPoolSize;
%                     'off' leaves the pool untouched, so run.run_fooof opens the
%                     profile default on its first parfor; a positive integer sets it
%                     explicitly. 'auto' only ever shrinks an existing pool, never grows
%                     one, so a pool you opened deliberately is left as you set it.
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
%   Figures
%   -------
%   backgroundfigs    Draw every figure off-screen: they are still created, printed and
%                     saved exactly as before, they just never appear on screen. Worth
%                     setting for a long batch, or on a machine you are using for
%                     something else while it runs - hundreds of figures popping up and
%                     stealing focus is the only thing this changes. They cannot be
%                     inspected interactively while off-screen. Default: false.
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

    %--- Parallel machines ---
    opts.uselocks (1,1) logical = true
    opts.lockpath         char  = ''
    opts.lockstalemin     (1,1) double {mustBePositive} = 360
    opts.lockheartbeatmin (1,1) double {mustBePositive} = 5

    %--- Parallel pool ---
    opts.poolworkers = 'auto'

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

    %--- Figures ---
    opts.backgroundfigs (1,1) logical = false


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
if isempty(opts.lockpath),     opts.lockpath     = fullfile(opts.figpath, '.locks'); end

%%% Draw off-screen, when asked. Set on the root's factory default rather than on each
%%% figure, so it covers every figure the plot functions and EEGLAB's topoplot create
%%% without any of them needing to know about it. print() renders an invisible figure
%%% exactly like a visible one, so the saved PNGs are identical either way; only the
%%% popping-up changes. The onCleanup puts the previous default back on every exit path,
%%% error and Ctrl-C included - this is global state, and leaving MATLAB in a mode where
%%% no figure ever shows again would be a nasty thing to hand back to the user.
if opts.backgroundfigs
    prevFigVisible  = get(groot, 'DefaultFigureVisible');
    restoreFigVisible = onCleanup(@() set(groot, 'DefaultFigureVisible', prevFigVisible));
    set(groot, 'DefaultFigureVisible', 'off');
    fprintf('Drawing figures off-screen (backgroundfigs=true).\n');
end

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

%%% Size the FOOOF pool before the batch rather than letting run.run_fooof open the
%%% profile default on its first parfor. That default is one worker per core whatever
%%% else is going on, and each worker is a separate MATLAB process competing with a
%%% client already holding two full-night recordings and eval_clean's interpolated
%%% copies of them. gedai.autoPoolSize scales that to the machine, so the same call is
%%% right on a 64 GB box and on a 2 TB one.
if ischar(opts.poolworkers) || isstring(opts.poolworkers)
    switch lower(string(opts.poolworkers))
        case "off",  nWorkers = [];
        case "auto", nWorkers = gedai.autoPoolSize();
        otherwise
            error('bidsfun_evalfigs:poolworkers', ...
                'poolworkers must be ''auto'', ''off'', or a positive integer.');
    end
else
    mustBeInteger(opts.poolworkers); mustBePositive(opts.poolworkers);
    nWorkers = double(opts.poolworkers);
end

if ~isempty(nWorkers)
    %%% Shrink an oversized pool, open one when there is none, and otherwise leave a
    %%% smaller existing pool alone - somebody sized it that way on purpose.
    try
        pool = gcp('nocreate');
        if isempty(pool)
            fprintf('Opening parallel pool with %d worker(s).\n', nWorkers);
            parpool('Processes', nWorkers);
        elseif pool.NumWorkers > nWorkers
            fprintf('Resizing parallel pool from %d to %d worker(s).\n', ...
                pool.NumWorkers, nWorkers);
            delete(pool);
            parpool('Processes', nWorkers);
        else
            fprintf('Reusing parallel pool (%d worker(s)).\n', pool.NumWorkers);
        end
    catch ME
        %%% No Parallel Computing Toolbox, or the pool refused to start. FOOOF still
        %%% runs, just serially - not a reason to abandon the batch.
        warning('bidsfun_evalfigs:noPool', ...
            'Could not size the parallel pool (%s); continuing without one.', ME.message);
    end
end

%%% Which figure means "finished". run.eval_clean writes its figures in a fixed order,
%%% so the last one the current toggles will produce is the only file whose presence
%%% proves the recording ran to completion. Gating on a fixed name instead - the first
%%% figure written - marks a recording done the moment it starts producing output, so
%%% anything that died partway through is skipped forever with an incomplete set.
sentinelSuffix = gedai.lastEvalFigure(opts);
if isempty(sentinelSuffix)
    fprintf(['No figure with a predictable name is enabled; every recording will be ' ...
             'processed regardless of existing output.\n']);
else
    fprintf('Resume check keys on _%s.png\n', sentinelSuffix);
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
    pngFile = fullfile(opts.figpath, ['desc-' opts.afterdesc], subDir, [fileID '_' sentinelSuffix '.png']);
    if ~opts.refresh && ~isempty(sentinelSuffix) && isfile(pngFile)
        fprintf('[skip] output exists: %s\n', pngFile)
        continue
    end

    %%% Claim this recording, so a second machine walking the same list moves on to the
    %%% next one instead of redoing this. The lock file is created with an atomic
    %%% create-if-absent, so two machines reaching this line together cannot both win.
    %%% lockGuard holds the claim: it releases on success, on the error caught below, and
    %%% on Ctrl-C or any error that unwinds out of this function, so a recording is never
    %%% left claimed by a run that is no longer working on it. Hence the explicit clear at
    %%% the end of the iteration - the claim must not outlive its recording.
    lockGuard = [];  %#ok<NASGU> release any claim still held from the previous iteration
    if opts.uselocks
        [lockGuard, acquired, holder] = claimFile( ...
            fullfile(opts.lockpath, [fileID '_desc-' opts.afterdesc '.lock']), ...
            'stalemin', opts.lockstalemin, 'heartbeatmin', opts.lockheartbeatmin); %#ok<ASGLU>
        if ~acquired
            fprintf('[skip] %s: claimed by %s (pid %d) since %s\n', ...
                fileID, holder.host, holder.pid, holder.started)
            continue
        end
        %%% Re-check now that the claim is ours: the other machine may have finished this
        %%% recording and dropped its lock between our skip check above and here.
        if ~opts.refresh && ~isempty(sentinelSuffix) && isfile(pngFile)
            fprintf('[skip] output exists: %s\n', pngFile)
            continue
        end
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

    catch ME
        fprintf('[ERROR] %s: %s\n', fileID, ME.message);
        failures{end+1} = struct('fileID', fileID, 'message', ME.message, 'report', ME.getReport()); %#ok<AGROW>
    end

    %%% Release this recording's figures and memory on every path, success or failure,
    %%% before the next one starts allocating. Both EEG structs run to several GB on a
    %%% full night and run.eval_clean holds its own interpolated copies while it works,
    %%% so carrying them across the boundary hands the next recording a peak it may not
    %%% fit inside. close all has to sit out here too: while it lived at the end of the
    %%% try, a recording that errored left every figure it had already drawn open, and a
    %%% run of consecutive failures piled them up - which is precisely what a machine
    %%% running out of memory produces.
    close all force
    clear EEGraw EEGafter

    %%% Drop the claim, whether the recording succeeded or failed. A failed recording has
    %%% to become available again - it is exactly the one another machine (or a later run
    %%% of this one) should be free to retry.
    clear lockGuard
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
