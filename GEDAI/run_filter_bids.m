function failures = run_filter_bids(BIDS, opts)
% RUN_FILTER_BIDS  Preprocess BIDS EEG files: import, resample, DC removal, bad
%   channel removal, Zapline.
%   Results are saved as EEGLAB .set files under
%   <BIDS root>/derivatives/prep-ged/<sub>/<ses>/. .set rather than BrainVision
%   because chanlocs and urchanlocs have to survive to the next stage: bad channels
%   are dropped here, and only urchanlocs says which ones to interpolate back.
%
% USAGE:
%   run_filter_bids(BIDS)
%   run_filter_bids(BIDS, subjectfilter={'sub-xxx'}, refresh=true)
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
%   badchannels     detect and remove bad channels after DC removal and before
%                   Zapline, while the line noise clean_channels keys on is still
%                   there. The mask is cached as
%                   <fileID>_desc-<desc>_badchans.mat next to the filtered file and
%                   read back by run_gedai_bids for the leadfield  (default true)
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
    error('run_filter_bids:noFiles', 'No matching EEG files found in BIDS layout.');
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
    EEG = eeg_import(eegFile);
    KeepTime = struct('EEGimport', toc(D));

    %%% Drop non-EEG channels
    %%% Done here rather than inside run.run_filter so the channel locations below line
    %%% up with the EEG channels (run.run_filter repeats it as a no-op).
    EEG = pop_select(EEG, 'nochannel', intersect(1:EEG.nbchan, opts.noteegchannels));

    %%% Read SFP file (dome-solved channel locations)
    %%% clean_channels needs coordinates, and .set output carries chanlocs/urchanlocs
    %%% forward so the channels removed below can be interpolated back downstream.
    if ~isempty(opts.sfppath)
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

    %%% Bad channel cache
    badchanFile = '';
    if opts.badchannels
        badchanFile = fullfile(outDir, [fileID '_desc-' opts.desc '_badchans.mat']);
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
        'refresh',            opts.refresh);

    if opts.plotResults
        nm = strrep(get(gcf, 'Name'), ' ', '_');
        print(gcf, fullfile(figDir, [fileID '_zapline_' nm '.png']), '-dpng', '-r100');
        pause(3); close(gcf);
    end

    %%% Bad channel figures: where the bad channels are, and when they went bad
    if opts.badchannels
        badchanFigDir = fullfile(opts.figpath, 'badchans', subDir);
        if ~exist(badchanFigDir, 'dir'), mkdir(badchanFigDir); end
        gedai.plotBadChannels(EEG.etc.badchans.corr, EEG.etc.badchans.znoise, ...
            EEG.urchanlocs, ...
            fullfile(badchanFigDir, [fileID '_desc-' opts.desc '_BadChannelTopoplot.png']), ...
            EEG.etc.badchans.flatprop);
        gedai.plotBadChannelTime(EEG.etc.badchans.corr, EEG.etc.badchans.mask, ...
            fullfile(badchanFigDir, [fileID '_desc-' opts.desc '_BadChannelTimecourse.png']), ...
            'windowseconds', EEG.etc.badchans.params.windowSeconds * ...
                             EEG.etc.badchans.params.windowStride, ...
            'title', fileID);
    end

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
        print(gcf, fullfile(figDir, [fileID '_zapline2_' nm '.png']), '-dpng', '-r100');
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

    %%% JSON sidecar: timings plus the parameters each step actually ran with.
    %%% The structs are built by run.run_filter next to the calls that consume them, so
    %%% what is recorded here cannot drift from what was applied.
    [~, baseName] = fileparts(outFile);
    prepParams = struct();
    if isfield(EEG.etc, 'filterparams'), prepParams = EEG.etc.filterparams; end
    prepParams.targetSampleRate = opts.targetsrate;
    prepParams.removeDC         = opts.removeDC;
    prepParams.zeroPatchSeconds = opts.zeropatchseconds;
    sidecarjson(KeepTime, fullfile(outDir, [baseName '.json']), ...
        struct('PreprocessingParameters', prepParams));

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
    fid = fopen(fullfile(opts.savepath, 'failed_files_zapline.json'), 'w');
    fprintf(fid, '%s', jsonencode([failures{:}], 'PrettyPrint', true));
    fclose(fid);
end
end
