function failures = bidsfun_hp_zap_cleanline(BIDS, opts)
% RUN_FILTER_BIDS  Preprocess BIDS EEG files: import, resample, DC removal, bad
%   channel removal, Zapline.
%   Results are saved as EEGLAB .set files under
%   <BIDS root>/derivatives/prep-ged/<sub>/<ses>/. .set rather than BrainVision
%   because chanlocs and urchanlocs have to survive to the next stage: bad channels
%   are dropped here, and only urchanlocs says which ones to interpolate back.
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
%   zeropatchseconds  cut out all-zero patches (amplifier crash padding) longer than
%                   this many seconds before filtering, restore them before saving;
%                   0 = skip                                    (default 5)
%   badchannels     remove the bad channels before Zapline, so they cannot influence
%                   its spatial filters                          (default true)
%   badchandesc     desc of the mask written by bidsfun_detect_badchans. Kept separate
%                   from desc so that re-filtering and re-detecting do not invalidate
%                   each other. If no mask exists this stage detects one itself, using
%                   the same code path                           (default 'badchan')
%   badchanavgref   average-reference the data for the detection only and undo it
%                   afterwards, so a single-electrode reference cannot make the ring
%                   of channels around it look bad              (default true)
%   badchanstride   evaluate every Nth window in the correlation criterion; the
%                   criterion is a proportion of windows, so this costs precision, not
%                   correctness. 1 restores the original behaviour  (default 2)
%   flatthreshold   peak-to-peak in uV below which a 5-s window counts as flat; a
%                   channel flat for more than half the recording is removed. A dead
%                   electrode passes both clean_channels criteria (0/0 = NaN, and NaN
%                   fails every comparison), so it needs its own test (default 0.5)
%   sfppath         path passed to the SFP resolver; clean_channels needs channel
%                   locations                                   (default BIDS root)
%   savefileext     '.set' (EEGLAB) or anything else for BrainVision (default '.set')
%   subjectfilter   cell array of subject ID strings; {} = all subjects
%   sessionfilter   cell array of session ID strings; {} = all sessions

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
    opts.noisefreqs                      = 50
    opts.adaptiveNremove  (1,1) logical  = true
    opts.fixedNremove     (1,1) double   = 1
    opts.chunkLength      (1,1) double   = 300
    opts.plotResults      (1,1) logical  = true
    opts.zeropatchseconds (1,1) double   = 5

    %--- Bad channels ---
    opts.badchannels      (1,1) logical  = true
    opts.badchandesc      char           = 'badchan'
    opts.badchanavgref    (1,1) logical  = true
    opts.badchanstride    (1,1) double   = 2
    opts.flatthreshold    (1,1) double   = 0.5
    opts.sfppath          char           = BIDS.pth

    %--- Subject filter ---
    opts.subjectfilter    cell            = {}
    opts.sessionfilter    cell            = {}
end

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
    %%% Done here rather than inside run.run_filter so the channel locations below line
    %%% up with the EEG channels (run.run_filter repeats it as a no-op).
    EEG = pop_select(EEG, 'nochannel', intersect(1:EEG.nbchan, opts.noteegchannels));

    %%% Channel locations
    %%% clean_channels needs coordinates, and .set output carries chanlocs/urchanlocs
    %%% forward so the channels removed below can be interpolated back downstream.
    EEG = gedai.assignChanlocs(EEG, BIDS, opts.sfppath, eegFile, p, fileID);

    %%% Bad channel mask, owned by bidsfun_detect_badchans and keyed to its own desc so
    %%% that re-filtering and re-detecting are independent. smartcache loads it when it
    %%% exists, so this stage only computes a mask when the detection stage was skipped.
    badchanFile = '';
    if opts.badchannels
        badchanFile = fullfile(outDir, [fileID '_desc-' opts.badchandesc '_badchans.mat']);
    end

    %%% Run filter pipeline
    [EEG, KeepTime] = run.run_filter(EEG, ...
        'noteegchannels', opts.noteegchannels, ...
        'targetsrate',    opts.targetsrate, ...
        'removeDC',       opts.removeDC, ...
        'zapline',        opts.zapline, ...
        'cleanline',      opts.cleanline, ...
        'KeepTime',        KeepTime, ...
        'noisefreqs',      opts.noisefreqs, ...
        'adaptiveNremove', opts.adaptiveNremove, ...
        'fixedNremove',    opts.fixedNremove, ...
        'chunkLength',     opts.chunkLength, ...
        'plotResults',     opts.plotResults, ...
        'zeropatchseconds',   opts.zeropatchseconds, ...
        'restorezeropatches', false, ...
        'badchannels',        opts.badchannels, ...
        'badchanavgref',      opts.badchanavgref, ...
        'badchanstride',      opts.badchanstride, ...
        'flatthreshold',      opts.flatthreshold, ...
        'badchanfile',        badchanFile, ...
        'refresh',            false);   % never re-detect here; that is the detection stage's job

    if opts.plotResults
        nm = strrep(get(gcf, 'Name'), ' ', '_');
        gedai.printFigure(gcf, fullfile(figDir, [fileID '_zapline_' nm '.png']));
        pause(3); close(gcf);
    end

    %%% The bad channel figures and channels.tsv belong to bidsfun_detect_badchans, which
    %%% owns that stage; drawing them again here would only duplicate them under a second
    %%% desc. The applied mask is still recorded in this stage's JSON sidecar.

    %%% Optional second Zapline pass
    if opts.zapline2
        D = tic; fprintf('\nZapline plus (pass 2) ...\n')
        [EEG.data, ~, ~] = clean_data_with_zapline_plus( ...
            double(EEG.data), EEG.srate, ...
            'noisefreqs',      opts.noisefreqs, ...
            'adaptiveNremove', opts.adaptiveNremove, ...
            'fixedNremove',    opts.fixedNremove, ...
            'chunkLength',     opts.chunkLength, ...
            'plotResults',     opts.plotResults);
        KeepTime.Zapline2 = toc(D);
        fprintf('ZapLine-plus pass 2: %.2f min\n', KeepTime.Zapline2 / 60);
    end
    if opts.plotResults & opts.zapline2
        nm = strrep(get(gcf, 'Name'), ' ', '_');
        gedai.printFigure(gcf, fullfile(figDir, [fileID '_zapline2_' nm '.png']));
        pause(3); close(gcf);
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
    %%% The structs are built by run.run_filter next to the calls that consume them, so
    %%% what is recorded here cannot drift from what was applied.
    [~, baseName] = fileparts(outFile);
    prepParams = struct();
    if isfield(EEG.etc, 'filterparams'), prepParams = EEG.etc.filterparams; end

    %%% Everything about bad channels goes to the badchans sidecar rather than this one,
    %%% so the mask, the figure and the settings that produced them are described in one
    %%% file. This stage's JSON is about filtering.
    bcParams = struct();
    if isfield(prepParams, 'BadChannels')
        bcParams   = prepParams.BadChannels;
        prepParams = rmfield(prepParams, 'BadChannels');
    end
    bcKeys   = intersect(fieldnames(KeepTime), ...
        {'FlatChannelDetection', 'BadChannelDetection'});
    bcTime   = struct();
    for k = 1:numel(bcKeys)
        bcTime.(bcKeys{k}) = KeepTime.(bcKeys{k});
    end
    KeepTime = rmfield(KeepTime, bcKeys);

    if opts.badchannels
        %%% Which mask the first round came from: this stage loads it rather than
        %%% detecting, so the sidecar has to name the file that decided it.
        bcParams.firstRoundDesc = opts.badchandesc;
        sidecarjson(bcTime, ...
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
        %%% Close whatever this iteration left open. Zapline's figure is created inside
        %%% run.run_filter but printed and closed here, so any failure in between - a
        %%% CleanLine error, a failed save - would otherwise orphan one figure per
        %%% recording, and a batch of nights ends with a screen full of them.
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
