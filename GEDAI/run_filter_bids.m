function run_filter_bids(BIDS, opts)
% RUN_FILTER_BIDS  Preprocess BIDS EEG files: import, resample, DC removal, Zapline.
%   Results are cached as BrainVision files under
%   <BIDS root>/derivatives/preprocessing/A_filtered/<sub>/<ses>/.
%
% USAGE:
%   run.run_filter_bids(BIDS)
%   run.run_filter_bids(BIDS, subjectfilter={'sub-xxx'}, refresh=true)
%
% INPUTS:
%   BIDS   — bids.layout object
%
% OPTIONAL NAME-VALUE:
%   savepath        output root directory
%                   (default <BIDS root>/derivatives/preprocessing/A_filtered)
%   refresh         force reprocessing even if cache file exists (default false)
%   tasklabel       BIDS task label(s) to query           (default {'Sleep','sleep'})
%   recordinglabel  BIDS recording label to query         (default '125Hz')
%   noteegchannels  channel indices to drop               (default 257:264)
%   targetsrate     resample target in Hz; 0 = skip       (default 125)
%   removeDC        apply DC-removal filter               (default true)
%   removeLN        apply Zapline-plus line-noise removal (default true)
%   subjectfilter   cell array of subject ID strings; {} = all subjects

arguments
    BIDS

    %--- Paths ---
    opts.savepath         char    = fullfile(BIDS.pth, 'derivatives', 'preprocessing')
    opts.bidsout (1,1)    logical = true
    opts.refresh (1,1)    logical = false
    opts.desc             char    = 'filt'

    %--- EEG ---
    opts.tasklabel                       = {'Sleep', 'sleep'}
    opts.recordinglabel   char           = '125Hz'
    opts.noteegchannels   (1,:) double   = 257:264
    opts.targetsrate      (1,1) double   = 125
    opts.removeDC         (1,1) logical  = true
    opts.removeLN         (1,1) logical  = true

    %--- Subject filter ---
    opts.subjectfilter    cell            = {}
end

%%% Query EEG files
filesEEG = bids.query(BIDS, 'data', 'extension', '.vhdr', ...
    'task', opts.tasklabel, 'recording', opts.recordinglabel);
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

    %%% Build output path
    if opts.bidsout
        outDir  = fullfile(opts.savepath, p.entities.sub, p.entities.ses);
        outFile = fullfile(outDir, [fileID '_desc-' opts.desc '_eeg.dat']);
    else
        outDir  = fullfile(opts.savepath);
        outFile = fullfile(outDir, 'A_filtered', [fileID '.mat']);
    end

    fprintf('Output → %s\n', outFile)

    %%% Skip if already cached and refresh not requested
    if ~opts.refresh && isfile(outFile)
        fprintf('[cached] skipping\n')
        continue
    end

    %%% Create output directory if needed
    outFileDir = fileparts(outFile);
    if ~exist(outFileDir, 'dir'), mkdir(outFileDir); end

    %%% Run filter pipeline
    [EEG, KeepTime] = run.run_filter(eegFile, ...
        'noteegchannels', opts.noteegchannels, ...
        'targetsrate',    opts.targetsrate, ...
        'removeDC',       opts.removeDC, ...
        'removeLN',       opts.removeLN);

    %%% Save EEG
    if opts.bidsout
        pop_writebva(EEG, outFile, 'DataOrientation', 'MULTIPLEXED');

        [~, baseName] = fileparts(outFile);
        sidecarjson(KeepTime, fullfile(outFileDir, [baseName '.json']));
    else
        save(outFile, 'EEG', 'KeepTime');
    end
end
end
