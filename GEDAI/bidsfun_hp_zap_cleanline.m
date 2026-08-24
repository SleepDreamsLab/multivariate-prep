function failures = bidsfun_hp_zap_cleanline(BIDS, opts)
% RUN_FILTER_BIDS  Preprocess BIDS EEG files: import, resample, DC removal, bad
%   channel removal, Zapline.
%   Results are saved as EEGLAB .set files under
%   <BIDS root>/derivatives/prep-ged/<sub>/<ses>/. .set rather than BrainVision
%   because chanlocs and urchanlocs have to survive to the next stage: bad channels
%   are dropped here, and only urchanlocs says which ones to interpolate back.
%
%   Line-noise removal (Zapline-plus, CleanLine) lives entirely in this file rather than
%   behind a shared filtering helper. Bad channels are NOT detected here: this stage only
%   loads the mask that bidsfun_detect_badchans already computed and cached. Running
%   detection twice, on data that has already been resampled/DC-filtered differently
%   between runs, is exactly the drift the two-stage split exists to avoid - so a missing
%   mask is a stop condition here, not a fallback path (see badchanmustexist below).
%
% USAGE:
%   bidsfun_hp_zap_cleanline(BIDS)
%   bidsfun_hp_zap_cleanline(BIDS, subjectfilter={'sub-xxx'}, refresh=true)
%
% INPUTS:
%   BIDS   — bids.layout object
%
% OPTIONAL NAME-VALUE:
%   savepath        output root directory
%                   (default <BIDS root>/derivatives/prep-ged)
%   refresh         force reprocessing even if output file exists (default false)
%   desc            BIDS desc entity for output filename         (default 'filt')
%   tasklabel       BIDS task label(s) to query                 (default {'Sleep','sleep'})
%   acqlabel        BIDS recording label to query               (default '125Hz')
%   noteegchannels  channel indices to drop                     (default 257:300)
%   targetsrate     resample target in Hz; 0 = skip             (default 125)
%   removeDC        apply DC-removal filter                     (default true)
%   zapline         apply Zapline-plus line-noise removal       (default true)
%   cleanline       apply CleanLine after Zapline               (default true)
%   zapline2        apply a second Zapline-plus pass after CleanLine (default false)
%   zapDetectionWinsize  window size in Hz for Zapline-plus' noise-peak detection,
%                   passed through as its detectionWinsize argument   (default 6)
%   notchfilt       apply an IIR notch filter at 50 and 100 Hz (filterbank's
%                   EEG_NotchFilt_IIR2) after Zapline/CleanLine, before the zero
%                   patches are restored                              (default false)
%   zeropatchseconds  cut out all-zero patches (amplifier crash padding) longer than
%                   this many seconds before filtering, restore them before saving;
%                   0 = skip                                    (default 5)
%   badchannels     remove the bad channels flagged by bidsfun_detect_badchans, before
%                   Zapline, so they cannot influence its spatial filters (default true)
%   badchandesc     desc of the mask written by bidsfun_detect_badchans. Kept separate
%                   from desc so that re-filtering and re-detecting do not invalidate
%                   each other                                   (default 'badchan')
%   badchanmustexist  when badchannels is true and no mask file exists yet: error if
%                   true, warn and proceed without removing any channel if false. Off by
%                   default it would be easy to silently filter a night nobody has
%                   screened yet                                 (default true)
%   sfppath         path passed to the SFP resolver; clean_channels needs channel
%                   locations                                   (default BIDS root)
%   savefileext     '.set' (EEGLAB) or anything else for BrainVision (default '.set')
%   subjectfilter   cell array of subject ID strings; {} = all subjects
%   sessionfilter   cell array of session ID strings; {} = all sessions
%
% Methods section:
%
% Continuous EEG (256 channels, 250 Hz) was high-pass filtered to remove DC offset and then
% cleaned of power-line artifacts using Zapline-plus (Klug & Kloosterman, 2022), an extension
% of the Zapline algorithm (de Cheveigné, 2020). Zapline separates the data into a
% line-frequency-dominated and a residual component using a comb filter, isolates the artifactual
% subspace via DSS applied to the line frequency and its harmonics, and removes it by spatial
% projection, thereby avoiding the spectral distortion introduced by notch filtering. The target
% frequency was fixed at 50 Hz; at a 250 Hz sampling rate the fundamental and its second harmonic
% (100 Hz) were removed jointly. To accommodate non-stationarity of the line artifact across the
% [8-h] recordings, data were processed in fixed 300-s chunks, with the individual noise peak
% re-estimated within each chunk (search window 50 ± 3 Hz) and the number of removed components
% determined adaptively per chunk by iterative outlier detection on the DSS component scores
% (minimum 1 component). The outlier-detection threshold was held fixed at sigma = 5, the upper
% bound of the algorithm's adaptive range and therefore its most conservative setting; the
% iterative threshold adaptation was disabled because it converged to this bound on every pilot
% recording, so a fixed threshold yields identical output while also equalising cleaning
% strength across recordings. Residual artifacts at 50 and 100 Hz were subsequently attenuated using
% a custom parallelised reimplementation CleanLine (Mullen, 2012; EEGLAB), leveraging
% multi-taper sinusoidal regression, .
% Within 4-s windows advanced in 2-s steps, the Thomson F-test for a deterministic sinusoid
% was evaluated at 50 and 100 Hz using 7 Slepian tapers (time–bandwidth product 4,
% resolution bandwidth 2 Hz) and an FFT zero-padded to 4096 points (0.061 Hz bin spacing).
% Sinusoids reaching p < 0.01 were fitted by least squares and subtracted; target frequencies
% not reaching significance in a given window and channel were left unmodified. Fitted
% components were crossfaded sigmoidally across the 2-s overlap between adjacent windows.
% Unlike the default CleanLine implementation, which subtracts an estimated sinusoid at each
% target frequency regardless of significance, subtraction here was conditional on the F-test.
%
% Bad channels, flagged by bidsfun_detect_badchans and removed here before Zapline-plus so
% that they cannot contribute to the estimation of its spatial filters (see that file for the
% detection criteria). The full montage was retained alongside the data so that removed channels
% could be restored by spherical-spline interpolation wherever whole-head output was required.

arguments
    BIDS

    %--- Paths ---
    opts.derivfolder      char    = 'prep-ged'
    opts.savepath         char    = ''
    opts.figpath          char    = ''
    opts.refresh (1,1)    logical = false
    opts.desc             char    = 'filt'
    opts.savefileext      char    = '.set'

    %--- EEG ---
    opts.tasklabel                       = {'Sleep', 'sleep'}
    opts.acqlabel   char                 = ''
    opts.noteegchannels   (1,:) double   = 257:300
    opts.targetsrate      (1,1) double   = 0
    opts.removeDC         (1,1) logical  = true
    opts.zapline          (1,1) logical  = true
    opts.cleanline        (1,1) logical  = true
    opts.zapline2         (1,1) logical  = false
    opts.notchfilt        (1,1) logical  = false
    opts.noisefreqs                      = 50
    opts.adaptiveNremove  (1,1) logical  = true
    opts.fixedNremove     (1,1) double   = 1
    opts.chunkLength      (1,1) double   = 0
    opts.plotResults      (1,1) logical  = true
    opts.zeropatchseconds (1,1) double   = 5
    opts.noiseCompDetectSigma         (1,1) double   = 3
    opts.adaptiveSigma    (1,1) logical  = true
    opts.zapDetectionWinsize (1,1) double = 6

    opts.cleanlineWinstep (1,1) double = 1
    opts.cleanlinePad (1,1) double = 0

    %--- Bad channels ---
    opts.badchannels      (1,1) logical  = true
    opts.badchandesc      char           = 'badchan'
    opts.badchanmustexist (1,1) logical  = true
    opts.sfppath          char           = BIDS.pth

    %--- Subject filter ---
    opts.subjectfilter    cell            = {}
    opts.sessionfilter    cell            = {}
end

fprintf('\n=== Running bidsfun_hp_zap_cleanline ===\n');

if isempty(opts.savepath), opts.savepath = fullfile(BIDS.pth, 'derivatives', opts.derivfolder); end
if isempty(opts.figpath),  opts.figpath  = fullfile(BIDS.pth, 'derivatives', opts.derivfolder, 'figures'); end

%%% Query EEG files
filesEEG = bids.query(BIDS, 'data', 'extension', '.vhdr', ...
    'task', opts.tasklabel, 'acq', opts.acqlabel);
if isempty(filesEEG)
    error('bidsfun_hp_zap_cleanline:noFiles', 'No matching EEG files found in BIDS layout.');
end

%%% Loop over EEG files
failures = {};
for ifile = 1:numel(filesEEG)
    eegFile = filesEEG{ifile};
    p       = bids.internal.parse_filename(eegFile);
    fileID  = strjoin(cellfun(@(k) [k '-' p.entities.(k)], fieldnames(p.entities), 'uni', 0), '_');

    %%% Subject filter
    if ~isempty(opts.subjectfilter) && ~contains(fileID, opts.subjectfilter)
        continue
    end

    %%% Session filter
    if ~isempty(opts.sessionfilter) && ~contains(fileID, opts.sessionfilter)
        continue
    end
    fprintf('\n=== %s ===\n', fileID)
    figsBefore = findall(0, 'Type', 'figure');
    try

    %%% Build output paths
    subDir   = fullfile(['sub-' p.entities.sub], ['ses-' p.entities.ses]);
    outDir   = fullfile(opts.savepath, subDir);
    outFile  = fullfile(outDir, [fileID '_desc-' opts.desc '_eeg' opts.savefileext]);
    figDir   = fullfile(opts.figpath, ['desc-' opts.desc], subDir);
    if ~exist(figDir, 'dir'), mkdir(figDir); end
    fprintf('Output → %s\n', outFile)

    %%% Skip if already processed and refresh not requested
    if ~opts.refresh && isfile(outFile)
        fprintf('[File already exists] skipping\n')
        continue
    end

    %%% Create output directory if needed
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    %%% Import EEG
    D = tic; fprintf('\nEEG import ...\n')
    EEG = fast_eeg_import(eegFile);
    KeepTime = struct('EEGimport', toc(D));

    %%% Drop non-EEG channels
    %%% Done here rather than inside run.run_prep so the channel locations below line
    %%% up with the EEG channels (run.run_prep repeats it as a no-op).
    EEG = pop_select(EEG, 'nochannel', intersect(1:EEG.nbchan, opts.noteegchannels));

    %%% Channel locations
    %%% clean_channels needs coordinates, and .set output carries chanlocs/urchanlocs
    %%% forward so the channels removed below can be interpolated back downstream.
    EEG = gedai.assignChanlocs(EEG, BIDS, opts.sfppath, eegFile, p, fileID);

    %%% Shared prep: resample, excise zero patches, DC removal. Identical to what
    %%% bidsfun_detect_badchans applied when it detected the mask this stage is about
    %%% to load, or the two stages would disagree about what data the mask describes.
    [EEG, KeepTime] = run.run_prep(EEG, ...
        'noteegchannels',   opts.noteegchannels, ...
        'targetsrate',      opts.targetsrate, ...
        'removeDC',         opts.removeDC, ...
        'zeropatchseconds', opts.zeropatchseconds, ...
        'KeepTime',         KeepTime);

    %%% Bad channels: load only. Detection is bidsfun_detect_badchans' job; this stage
    %%% never redetects, so a stale or missing mask surfaces immediately rather than
    %%% silently filtering a night nobody has looked at yet.
    removedMask = false(EEG.nbchan, 1);
    badchanFile = '';
    if opts.badchannels
        badchanFile = fullfile(outDir, [fileID '_desc-' opts.badchandesc '_badchans.mat']);
        if isfile(badchanFile)
            S = load(badchanFile, 'removed_channels');
            removedMask = S.removed_channels(:);
        else
            msg = sprintf(['No bad-channel mask found at %s. Run bidsfun_detect_badchans ' ...
                'first, or set badchannels=false to skip removal.'], badchanFile);
            if opts.badchanmustexist
                error('bidsfun_hp_zap_cleanline:noBadChanFile', '%s', msg);
            else
                warning('bidsfun_hp_zap_cleanline:noBadChanFile', '%s', msg);
            end
        end
    end

    if opts.badchannels
        %%% Keep a record of the full montage before anything is dropped; urchanlocs is
        %%% how bidsfun_gedai later works out which channels to interpolate back.
        if ~isfield(EEG, 'urchanlocs') || isempty(EEG.urchanlocs)
            EEG.urchanlocs = EEG.chanlocs;
            for iCh = 1:numel(EEG.chanlocs)
                EEG.chanlocs(iCh).urchan = iCh;
            end
        end
        fprintf('Removing %d/%d channels flagged as bad.\n', nnz(removedMask), EEG.nbchan)
        EEG = pop_select(EEG, 'nochannel', find(removedMask));
        EEG.etc.filterparams.BadChannels = struct( ...
            'maskFile', badchanFile, 'nRemoved', nnz(removedMask));
    end

    %%% Zapline
    if opts.zapline
        D = tic; fprintf('\nZapline plus ...\n')
        %%% Cast to double up front: the old single-precision buffer (~7 GB at 256 ch
        %%% x 8 h) has no references left once this assignment completes and is freed
        %%% before the call, rather than sitting resident alongside the double copy
        %%% for the whole (multi-minute) duration of clean_data_with_zapline_plus.
        EEG.data = double(EEG.data);
        %%% Transpose to [nSamples x nChannels] here rather than letting
        %%% clean_data_with_zapline_plus transpose it internally. MATLAB cannot release
        %%% this EEG.data field's pre-call value until the assignment below completes,
        %%% so an internal transpose would sit resident as a second full copy of the
        %%% recording (~15 GB at 256 ch x 8 h) for the entire multi-minute call.
        %%% Transposing here instead pays that doubled-memory peak only briefly, before
        %%% the call starts, since clean_data_with_zapline_plus auto-detects that data
        %%% already has more rows than columns and skips its own transpose.
        EEG.data = EEG.data.';
        [EEG.data, zaplineConfig, analyticsResults] = clean_data_with_zapline_plus( ...
            EEG.data, EEG.srate, ...
            'noisefreqs',      opts.noisefreqs, ...
            'adaptiveNremove', opts.adaptiveNremove, ...
            'fixedNremove',    opts.fixedNremove, ...
            'chunkLength',     opts.chunkLength, ...
            'noiseCompDetectSigma', opts.noiseCompDetectSigma, ...
            'adaptiveSigma',   opts.adaptiveSigma, ...
            'detectionWinsize', opts.zapDetectionWinsize, ...
            'plotResults',     opts.plotResults);
        EEG.data = EEG.data.';
        EEG.etc.zapline.config    = zaplineConfig;
        EEG.etc.zapline.analytics = analyticsResults;
        EEG.etc.filterparams.Zapline = struct( ...
            'noisefreqs',           opts.noisefreqs, ...
            'adaptiveNremove',      opts.adaptiveNremove, ...
            'fixedNremove',         opts.fixedNremove, ...
            'chunkLengthSeconds',   opts.chunkLength, ...
            'noiseCompDetectSigma', opts.noiseCompDetectSigma, ...
            'adaptiveSigma',        opts.adaptiveSigma, ...
            'detectionWinsize',     opts.zapDetectionWinsize);
        KeepTime.Zapline = toc(D);
        fprintf('ZapLine-plus: %.2f min\n', KeepTime.Zapline / 60);

        if opts.plotResults
            nm = strrep(get(gcf, 'Name'), ' ', '_');
            gedai.printFigure(gcf, fullfile(figDir, [fileID '_zapline_' nm '.png']));
            pause(3); close(gcf);
        end
    end

    %%% Cleanline
    if opts.cleanline
        D = tic; fprintf('\nClean line ...\n')
        clp = struct('linefreqs', [opts.noisefreqs opts.noisefreqs*2], ...
            'winsize', 4, 'winstep', opts.cleanlineWinstep, 'sigtest', true, 'pad', opts.cleanlinePad);
        EEG = cleanline_fast(EEG, 'linefreqs', clp.linefreqs, ...
            'winsize', clp.winsize, 'winstep', clp.winstep, ...
            'sigtest', clp.sigtest, 'pad', clp.pad);
        EEG.etc.filterparams.CleanLine = clp;
        KeepTime.Cleanline = toc(D);
    end

    %%% The bad channel figures and channels.tsv belong to bidsfun_detect_badchans, which
    %%% owns that stage; drawing them again here would only duplicate them under a second
    %%% desc. The applied mask is still recorded in this stage's JSON sidecar.

    %%% Optional second Zapline pass
    if opts.zapline2
        D = tic; fprintf('\nZapline plus (pass 2) ...\n')
        EEG.data = double(EEG.data);
        %%% See the memory comment on the first Zapline pass above: pre-transposing
        %%% here avoids clean_data_with_zapline_plus holding a second full copy of the
        %%% recording resident for the whole call.
        EEG.data = EEG.data.';
        [EEG.data, ~, ~] = clean_data_with_zapline_plus( ...
            EEG.data, EEG.srate, ...
            'noisefreqs',      opts.noisefreqs, ...
            'adaptiveNremove', opts.adaptiveNremove, ...
            'fixedNremove',    opts.fixedNremove, ...
            'chunkLength',     opts.chunkLength, ...
            'detectionWinsize', opts.zapDetectionWinsize, ...
            'plotResults',     opts.plotResults);
        EEG.data = EEG.data.';
        KeepTime.Zapline2 = toc(D);
        fprintf('ZapLine-plus pass 2: %.2f min\n', KeepTime.Zapline2 / 60);

        if opts.plotResults
            nm = strrep(get(gcf, 'Name'), ' ', '_');
            gedai.printFigure(gcf, fullfile(figDir, [fileID '_zapline2_' nm '.png']));
            pause(3); close(gcf);
        end
    end

    %%% Optional notch filter at 50/100 Hz, after Zapline/CleanLine
    if opts.notchfilt
        D = tic; fprintf('\nNotch filter (50/100 Hz) ...\n')
        NotchFilters = filterbank(EEG.srate, 'EEG_NotchFilt_IIR2');
        for iFilt = 1:numel(NotchFilters)
            EEG.data = filtfilt(NotchFilters{iFilt}, double(EEG.data)')';
        end
        KeepTime.NotchFilter = toc(D);
    end

    %%% Put the all-zero patches back, so the saved file keeps its original length
    EEG = run.restore_zero_patches(EEG);

    %%% Save output
    %%% .set by default: BrainVision stores neither chanlocs nor urchanlocs, and both
    %%% are needed downstream to interpolate the channels removed as bad.
    EEG.data = single(EEG.data);
    if strcmpi(opts.savefileext, '.set')
        [~, outName, outExt] = fileparts(outFile);
        pop_saveset(EEG, 'filename', [outName outExt], 'filepath', outDir);
    else
        pop_writebva(EEG, outFile, 'DataOrientation', 'MULTIPLEXED');
    end

    %%% JSON sidecars: timings plus the parameters each step actually ran with.
    [~, baseName] = fileparts(outFile);
    prepParams = struct();
    if isfield(EEG.etc, 'filterparams'), prepParams = EEG.etc.filterparams; end

    %%% Everything about bad channels goes to the badchans sidecar rather than this one,
    %%% so which mask was applied is described in its own file. This stage's JSON is
    %%% about filtering.
    bcParams = struct();
    if isfield(prepParams, 'BadChannels')
        bcParams   = prepParams.BadChannels;
        prepParams = rmfield(prepParams, 'BadChannels');
    end

    if opts.badchannels
        bcParams.firstRoundDesc = opts.badchandesc;
        sidecarjson(struct(), ...
            fullfile(outDir, [fileID '_desc-' opts.desc '_badchans.json']), ...
            struct('BadChannelParameters', bcParams));
    end

    prepParams.targetSampleRate = opts.targetsrate;
    prepParams.removeDC         = opts.removeDC;
    prepParams.zeroPatchSeconds = opts.zeropatchseconds;
    sidecarjson(KeepTime, fullfile(outDir, [baseName '.json']), ...
        struct('PreprocessingParameters', prepParams));

    catch ME
        fprintf('[ERROR] %s: %s\n', fileID, ME.message);
        %%% Close whatever this iteration left open. Zapline's figure is printed and
        %%% closed right after it is created, but any failure downstream - a CleanLine
        %%% error, a failed save - would otherwise orphan figures from later steps, and
        %%% a batch of nights ends with a screen full of them.
        figsNow = findall(0, 'Type', 'figure');
        close(figsNow(~ismember(figsNow, figsBefore)));
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
    fid = fopen(fullfile(opts.savepath, 'failed_files_zapline.json'), 'w');
    fprintf(fid, '%s', jsonencode([failures{:}], 'PrettyPrint', true));
    fclose(fid);
end
end
