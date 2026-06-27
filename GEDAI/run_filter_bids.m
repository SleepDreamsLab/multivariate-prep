function run_filter_bids(BIDS, opts)
% RUN_FILTER_BIDS  Preprocess BIDS EEG files: import, resample, DC removal, Zapline.
%   Results are saved as BrainVision files under
%   <BIDS root>/derivatives/preprocessing/<sub>/<ses>/.
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
%                   (default <BIDS root>/derivatives/preprocessing)
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
%   subjectfilter   cell array of subject ID strings; {} = all subjects

arguments
    BIDS

    %--- Paths ---
    opts.savepath         char    = fullfile(BIDS.pth, 'derivatives', 'preprocessing')
    opts.refresh (1,1)    logical = false
    opts.desc             char    = 'filt'

    %--- EEG ---
    opts.tasklabel                       = {'Sleep', 'sleep'}
    opts.acqlabel   char           = '125Hz'
    opts.noteegchannels   (1,:) double   = 257:264
    opts.targetsrate      (1,1) double   = 125
    opts.removeDC         (1,1) logical  = true
    opts.zapline          (1,1) logical  = true
    opts.cleanline        (1,1) logical  = false
    opts.zapline2         (1,1) logical  = true

    %--- Subject filter ---
    opts.subjectfilter    cell            = {}
end

%%% Query EEG files
filesEEG = bids.query(BIDS, 'data', 'extension', '.vhdr', ...
    'task', opts.tasklabel, 'acq', opts.acqlabel);
if isempty(filesEEG)
    error('run_filter_bids:noFiles', 'No matching EEG files found in BIDS layout.');
end

%%% Loop over EEG files
for ifile = 1:numel(filesEEG)
    eegFile = filesEEG{ifile};
    p       = bids.internal.parse_filename(eegFile);
    fileID  = strjoin(cellfun(@(k) [k '-' p.entities.(k)], fieldnames(p.entities), 'uni', 0), '_');

    %%% Subject filter
    if ~isempty(opts.subjectfilter) && ~contains(fileID, opts.subjectfilter)
        continue
    end
    fprintf('\n=== %s ===\n', fileID)

    %%% Build output paths
    outDir   = fullfile(opts.savepath, ['sub-' p.entities.sub], ['ses-' p.entities.ses]);
    outFile  = fullfile(outDir, [fileID '_desc-' opts.desc '_eeg.dat']);
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
        'KeepTime',       KeepTime);

    %%% Optional second Zapline pass
    if opts.zapline2
        D = tic; fprintf('\nZapline plus (pass 2) ...\n')
        [EEG.data, ~, ~] = clean_data_with_zapline_plus( ...
            double(EEG.data), EEG.srate, ...
            'noisefreqs', 'line', ...
            'maxfreq', floor(EEG.srate/2) - 5, ...
            'plotResults', 0);
        KeepTime.Zapline2 = toc(D);
        fprintf('ZapLine-plus pass 2: %.2f min\n', KeepTime.Zapline2 / 60);
    end

    %%% Save BrainVision output
    EEG.data = single(EEG.data);
    pop_writebva(EEG, outFile, 'DataOrientation', 'MULTIPLEXED');

    %%% JSON timing sidecar
    [~, baseName] = fileparts(outFile);
    sidecarjson(KeepTime, fullfile(outDir, [baseName '.json']));
end
end
