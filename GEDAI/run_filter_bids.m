function failures = run_filter_bids(BIDS, opts)
% RUN_FILTER_BIDS  Preprocess BIDS EEG files: import, resample, DC removal, Zapline.
%   Results are saved as BrainVision files under
%   <BIDS root>/derivatives/prep-ged/<sub>/<ses>/.
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
%   noteegchannels  channel indices to drop                     (default 257:264)
%   targetsrate     resample target in Hz; 0 = skip             (default 125)
%   removeDC        apply DC-removal filter                     (default true)
%   zapline         apply Zapline-plus line-noise removal       (default true)
%   cleanline       apply CleanLine after Zapline               (default true)
%   zapline2        apply a second Zapline-plus pass after CleanLine (default false)
%   zeropatchseconds  cut out all-zero patches (amplifier crash padding) longer than
%                   this many seconds before filtering, restore them before saving;
%                   0 = skip                                    (default 5)
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

    %--- EEG ---
    opts.tasklabel                       = {'Sleep', 'sleep'}
    opts.acqlabel   char                 = ''
    opts.noteegchannels   (1,:) double   = 257:299
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
    outFile  = fullfile(outDir, [fileID '_desc-' opts.desc '_eeg.dat']);
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
        'restorezeropatches', false);
    if opts.plotResults
        nm = strrep(get(gcf, 'Name'), ' ', '_');
        print(gcf, fullfile(figDir, [fileID '_zapline_' nm '.png']), '-dpng', '-r150');
        pause(3); close(gcf);
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
        print(gcf, fullfile(figDir, [fileID '_zapline2_' nm '.png']), '-dpng', '-r150');
        pause(3); close(gcf);
    end

    %%% Put the all-zero patches back, so the saved file keeps its original length
    EEG = run.restore_zero_patches(EEG);

    %%% Save BrainVision output
    EEG.data = single(EEG.data);
    pop_writebva(EEG, outFile, 'DataOrientation', 'MULTIPLEXED');

    %%% JSON timing sidecar
    [~, baseName] = fileparts(outFile);
    sidecarjson(KeepTime, fullfile(outDir, [baseName '.json']));

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
