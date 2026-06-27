function eval_filt_bids(BIDS, opts)
% EVAL_FILT_BIDS  Evaluate filtering quality by comparing raw and filtered EEG.
%
%   eval_filt_bids(BIDS)
%   eval_filt_bids(BIDS, Name, Value, ...)
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
%                     Default: <BIDS root>/derivatives/preprocessing
%   filtdesc          desc label for the filtered filename. Default: 'filt'
%   scoringpath       Directory containing sleep-scoring files.
%                     Default: <BIDS root>/derivatives/scoring/scores/Manual_Checked
%   sfppath           Path passed to the SFP resolver. Default: <BIDS root>
%
%   Output paths
%   ------------
%   figpath           Root for evaluation figures.
%                     Default: <BIDS root>/derivatives/preprocessing/figures
%   refresh           Force re-run even if cached files exist. Default: false.
%
%   EEG
%   ---
%   tasklabel         BIDS task label(s) to query. Default: {'Sleep','sleep'}.
%   recordinglabel    BIDS acq label to query. Default: '125Hz'.
%   noteegchannels    Channel indices to drop from raw EEG. Default: 257:264.
%   net               EEG net identifier passed to chans1020. Default: 'EGI256'.
%
%   Subject filter
%   --------------
%   subjectfilter     Cell array of subject ID strings; {} = all subjects.
%
%   Epoch
%   -----
%   epochlength       Epoch duration in seconds. Default: 30.

arguments
    BIDS

    %--- Input paths ---
    opts.filteredpath     char = fullfile(BIDS.pth, 'derivatives', 'preprocessing')
    opts.filtdesc         char = 'filt'
    opts.scoringpath      char = fullfile(BIDS.pth, 'derivatives', 'scoring', 'scores', 'Manual_Checked')
    opts.sfppath          char = BIDS.pth

    %--- Output paths ---
    opts.figpath          char = fullfile(BIDS.pth, 'derivatives', 'preprocessing', 'figures')
    opts.refresh (1,1)    logical = false

    %--- EEG ---
    opts.tasklabel                      = {'Sleep', 'sleep'}
    opts.recordinglabel   char         = '125Hz'
    opts.noteegchannels   (1,:) double = 257:264
    opts.net              char         = 'EGI256'

    %--- Subject filter ---
    opts.subjectfilter    cell          = {}

    %--- Epoch ---
    opts.epochlength (1,1) double = 30
end

%%% Query raw EEG files from BIDS
filesEEG = bids.query(BIDS, 'data', 'extension', '.vhdr', ...
    'task', opts.tasklabel, 'acq', opts.recordinglabel);
if isempty(filesEEG)
    error('eval_filt_bids:noFiles', 'No matching EEG files found in BIDS layout.');
end

%%% Scoring files
scoringfiles = gedai.collectScoringFiles(opts.scoringpath);

%%% Loop over EEG files
for ifile = 1:numel(filesEEG)
    rawFile = filesEEG{ifile};
    p       = bids.internal.parse_filename(rawFile);
    fileID  = strjoin(cellfun(@(k) [k '-' p.entities.(k)], fieldnames(p.entities), 'uni', 0), '_');
    subDir  = fullfile(['sub-' p.entities.sub], ['ses-' p.entities.ses]);

    %%% Subject filter
    if ~isempty(opts.subjectfilter) && ~contains(fileID, opts.subjectfilter)
        continue
    end
    fprintf('\n=== %s ===\n', fileID)

    %%% Resolve filtered file
    filtFile = fullfile(opts.filteredpath, subDir, [fileID '_desc-' opts.filtdesc '_eeg.vhdr']);
    if ~isfile(filtFile)
        fprintf('[skip] filtered file not found: %s\n', filtFile)
        continue
    end
    fprintf('Raw    → %s\nFilt   → %s\n', rawFile, filtFile)

    %%% Find and load scoring
    scoringFile = gedai.matchScoringFile(p.entities, scoringfiles);
    if isempty(scoringFile)
        error('eval_filt_bids:noScoring', 'No scoring file matched for %s.', fileID);
    end
    fprintf('Scoring → %s\n', scoringFile)
    scoringDigits = scoreloader(scoringFile);

    %%% Load SFP
    sfpFile      = gedai.matchSfpFile(opts.sfppath, p.entities.sub, p.entities.ses);
    fprintf('SFP     → %s\n', sfpFile)
    chanlocs_reg = register_fiducials(readlocs(sfpFile));

    %%% Import raw EEG
    fprintf('Importing raw EEG ...\n')
    EEGraw = eeg_import(rawFile);
    EEGraw = pop_select(EEGraw, 'nochannel', intersect(1:EEGraw.nbchan, opts.noteegchannels));
    EEGraw.chanlocs   = chanlocs_reg(1:EEGraw.nbchan);
    EEGraw.urchanlocs = EEGraw.chanlocs;

    %%% Import filtered EEG
    fprintf('Importing filtered EEG ...\n')
    EEGfilt = eeg_import(filtFile);
    EEGfilt.chanlocs   = chanlocs_reg(1:EEGfilt.nbchan);
    EEGfilt.urchanlocs = EEGfilt.chanlocs;

%     %%% Trim scoring to match filtered EEG epoch count
%     nEpochs = floor(EEGfilt.pnts / (opts.epochlength * EEGfilt.srate));
%     while numel(scoringDigits) > nEpochs; scoringDigits(end) = []; end
%     scoringDigits = gedai.killN1(scoringDigits);

    %%% Evaluate
    figDir = fullfile(opts.figpath, ['desc-' opts.filtdesc], subDir);
    if ~exist(figDir, 'dir'), mkdir(figDir); end

    run.eval_clean(EEGraw, EEGfilt, scoringDigits, ...
        'EpochLength', opts.epochlength, 'WelchWindow', 4, ...
        'SavePath', fullfile(figDir, fileID), ...
        'refresh', opts.refresh, 'net', opts.net);
    close all;
end
end
