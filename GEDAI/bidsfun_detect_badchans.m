function failures = bidsfun_detect_badchans(BIDS, opts)
% BIDSFUN_DETECT_BADCHANS  Detect flat and bad channels, before any filtering.
%   Writes, per recording, under <BIDS root>/derivatives/prep-ged/<sub>/<ses>/:
%     <fileID>_desc-<desc>_badchans.mat   mask + per-window correlations + znoise
%     <fileID>_desc-<desc>_channels.tsv   BIDS channel table with a status column
%     <fileID>_desc-<desc>_badchans.json  timings and the parameters used
%   plus topoplot and timecourse figures under <figpath>/badchans/<sub>/<ses>/.
%
%   This is its own stage, ahead of bidsfun_hp_zap_cleanline, for two reasons. Detection is
%   a decision made once per recording while filtering may be re-run repeatedly, and
%   keying the mask to its own desc keeps the two from invalidating each other. It is
%   also the step whose output you want to look at - and possibly overrule - before
%   committing to the expensive stages downstream.
%
%   Bad-channel detection lives entirely in this file rather than behind a shared
%   filtering helper, so it can never be pulled out of step with line-noise removal
%   (bidsfun_hp_zap_cleanline, which loads the .mat this stage writes). Only the generic
%   prep both stages need identically - drop channels, resample, excise zero patches,
%   DC removal - is shared, via run.run_prep.
%
% USAGE:
%   bidsfun_detect_badchans(BIDS)
%   bidsfun_detect_badchans(BIDS, subjectfilter={'sub-xxx'}, refresh=true)
%
% OPTIONAL NAME-VALUE:
%   desc            BIDS desc entity for the outputs             (default 'badchan')
%   refresh         redetect even if the .mat exists             (default false)
%   targetsrate     resample target in Hz; 0 = skip. Match bidsfun_hp_zap_cleanline, or
%                   detection and filtering see different data   (default 0)
%   removeDC        apply the DC-removal filter first. Likewise  (default true)
%   zeropatchseconds  cut all-zero patches longer than this before detecting (default 5)
%   badchanavgref       average-reference for detection                     (default true)
%   badchanstride       evaluate every Nth window in the correlation criterion (default 2)
%   badchancorrthreshold  reconstruction correlation below which a window fails the
%                       criterion                                          (default 0.8)
%   badchanmaxbrokentime  proportion of windows failing the correlation criterion
%                       above which a channel is removed                  (default 0.3)
%   badchannoisethreshold  robust z-score of the line-noise ratio above which a
%                       channel fails the noise criterion                 (default 4)
%   badchanwindowseconds  window length, in seconds, for both criteria    (default 5)
%   badchannumsamples   number of RANSAC subsets sampled for the correlation criterion
%                                                                          (default 50)
%   badchansubsetfraction  fraction of channels in each RANSAC subset     (default 0.25)
%   flatthreshold       peak-to-peak uV below which a window is flat      (default 0.5)
%   flatmaxbrokentime   proportion of windows failing the flat criterion above which
%                       a channel is removed                              (default 0.4)
%   noteegchannels      channel indices to drop                           (default 257:300)
%   ramsaver            filter one channel at a time instead of the whole data matrix
%                       at once (DC removal, the pre-RANSAC HP filter, and the noise
%                       criterion's >50Hz removal), trading speed for a lower peak
%                       memory footprint                                  (default false)
%   sfppath         path passed to the SFP resolver              (default BIDS root)
%   subjectfilter   cell array of subject ID strings; {} = all subjects
%   sessionfilter   cell array of session ID strings; {} = all sessions
%
% Bad channels:
%
% Bad channels were identified and removed before line-noise correction. For detection only, the
% recording was re-referenced to the common average; the data carried forward retained the original
% vertex (Cz) reference. This is necessary because under a single-electrode reference the electrodes
% immediately surrounding the reference carry a near-zero potential difference, and both criteria
% below are ratios in which that difference forms the denominator, so those channels are
% systematically misclassified (here, the median reconstruction correlation of the ring of
% electrodes adjacent to the reference fell from 0.99 to 0.73 when detection was run on
% vertex-referenced data, while every other channel was unaffected). Three criteria were evaluated
% over consecutive 5-s windows, and a channel was rejected when either the flat or the correlation
% criterion held for more than 30% of the recording. (i) Flat: the signal varied by no more than 0.5 uV peak-to-peak within the
% window. This criterion is required because a dead electrode is invisible to the two that follow -
% both are ratios that evaluate to 0/0 for a constant signal, and the resulting NaN fails every
% threshold comparison. (ii) Line noise: the ratio of the >45 Hz to the <45 Hz robust amplitude
% (median absolute deviation), the two bands separated by a 100th-order least-squares FIR lowpass
% with a 45-50 Hz transition, expressed as a robust z-score across channels and thresholded at
% z = 4. Because this criterion keys on power-line contamination, detection precedes Zapline-plus,
% while the 50 and 100 Hz components are still present. (iii) Reconstruction correlation: 25 random
% subsets of 25% of the channels were each used to predict all channels by spherical-spline
% interpolation (Perrin et al., 1989); the element-wise median of the 25 predictions was taken as
% the consensus estimate (RANSAC; Fischler & Bolles, 1981), and the Pearson correlation between each
% channel's <45 Hz signal and its consensus prediction computed per window, with a threshold of
% r = 0.8. Windows were evaluated at a stride of two, leaving the proportion estimated from
% approximately 3,300 windows per recording. Channels flagged by any criterion were removed before
% Zapline-plus so that they could not contribute to the estimation of its spatial filters; the full
% montage was retained alongside the data so that removed channels could be restored by
% spherical-spline interpolation wherever whole-head output was required. Detection used
% clean_channels from the clean_rawdata plugin for EEGLAB (Kothe & Makeig, 2013; Delorme & Makeig,
% 2004), locally modified to return the per-window correlations and to accept a window stride; the
% flat-line criterion was implemented separately.

arguments
    BIDS

    %--- Paths ---
    opts.derivfolder      char    = 'prep-ged'
    opts.savepath         char    = ''
    opts.figpath           char    = ''
    opts.refresh (1,1)    logical = false
    opts.desc              char    = 'badchan'

    %--- EEG ---
    opts.tasklabel                       = {'Sleep', 'sleep'}
    opts.acqlabel   char                 = ''
    opts.noteegchannels   (1,:) double   = 257:300
    opts.targetsrate      (1,1) double   = 0
    opts.removeDC         (1,1) logical  = true
    opts.zeropatchseconds (1,1) double   = 5
    opts.ramsaver         (1,1) logical  = false

    %--- Bad channels ---
    opts.badchanavgref        (1,1) logical  = true
    opts.badchanstride        (1,1) double   = 1
    opts.badchancorrthreshold (1,1) double   = 0.8
    opts.badchanmaxbrokentime (1,1) double   = 0.4
    opts.badchannoisethreshold (1,1) double  = 4
    opts.badchanwindowseconds (1,1) double   = 5
    opts.badchannumsamples    (1,1) double   = 50
    opts.badchansubsetfraction (1,1) double  = 0.25
    opts.badchanminbrokentime (1,1) double   = .10
    opts.flatthreshold        (1,1) double   = 0.5
    opts.flatmaxbrokentime    (1,1) double   = 0.4
    opts.sfppath               char           = BIDS.pth

    %--- Subject filter ---
    opts.subjectfilter    cell            = {}
    opts.sessionfilter    cell            = {}
end

fprintf('\n=== Running bidsfun_detect_badchans ===\n');

if isempty(opts.savepath), opts.savepath = fullfile(BIDS.pth, 'derivatives', opts.derivfolder); end
if isempty(opts.figpath),  opts.figpath  = fullfile(BIDS.pth, 'derivatives', opts.derivfolder, 'figures'); end

%%% Query EEG files
filesEEG = bids.query(BIDS, 'data', 'extension', '.vhdr', ...
    'task', opts.tasklabel, 'acq', opts.acqlabel);
if isempty(filesEEG)
    error('bidsfun_detect_badchans:noFiles', 'No matching EEG files found in BIDS layout.');
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
    subDir      = fullfile(['sub-' p.entities.sub], ['ses-' p.entities.ses]);
    outDir      = fullfile(opts.savepath, subDir);
    badchanFile = fullfile(outDir, [fileID '_desc-' opts.desc '_badchans.mat']);
    figDir      = fullfile(opts.figpath, 'badchans', subDir);
    fprintf('Output → %s\n', badchanFile)

    %%% Skip if already detected and refresh not requested
    if ~opts.refresh && isfile(badchanFile)
        fprintf('[File already exists] skipping\n')
        continue
    end

    if ~exist(outDir, 'dir'), mkdir(outDir); end
    if ~exist(figDir, 'dir'), mkdir(figDir); end

    %%% Import EEG
    D = tic; fprintf('\nEEG import ...\n')
    EEG = fast_eeg_import(eegFile);
    KeepTime = struct('EEGimport', toc(D));

    %%% Drop non-EEG channels, so the channel locations below line up
    EEG = pop_select(EEG, 'nochannel', intersect(1:EEG.nbchan, opts.noteegchannels));

    %%% Channel locations; clean_channels needs coordinates
    EEG = gedai.assignChanlocs(EEG, BIDS, opts.sfppath, eegFile, p, fileID);

    if isempty(EEG.chanlocs) || ~isfield(EEG.chanlocs, 'X') || isempty([EEG.chanlocs.X])
        error('bidsfun_detect_badchans:noChanlocs', ...
            'Bad channel detection needs EEG.chanlocs with X/Y/Z coordinates.')
    end

    %%% Shared prep: resample, excise zero patches, DC removal. Identical to what
    %%% bidsfun_hp_zap_cleanline applies later, so the data detection sees here is
    %%% exactly what filtering will see.
    [EEG, KeepTime] = run.run_prep(EEG, ...
        'noteegchannels',     opts.noteegchannels, ...
        'targetsrate',        opts.targetsrate, ...
        'removeDC',           opts.removeDC, ...
        'zeropatchseconds',   opts.zeropatchseconds, ...
        'ramsaver',           opts.ramsaver, ...
        'KeepTime',           KeepTime);

    %%% Keep a record of the full montage before anything is dropped
    if ~isfield(EEG, 'urchanlocs') || isempty(EEG.urchanlocs)
        EEG.urchanlocs = EEG.chanlocs;
        for iCh = 1:numel(EEG.chanlocs)
            EEG.chanlocs(iCh).urchan = iCh;
        end
    end

    %%% Parameters for both criteria, kept in one place: this is what gets written to
    %%% the JSON sidecar below, so what is recorded is necessarily what was run.
    bcp = struct( ...
        'corrThreshold',       opts.badchancorrthreshold, ...
        'noiseThreshold',      opts.badchannoisethreshold, ...
        'windowSeconds',       opts.badchanwindowseconds, ...
        'maxBrokenTime',       opts.badchanmaxbrokentime, ...
        'numSamples',          opts.badchannumsamples, ...
        'subsetSizeFraction',  opts.badchansubsetfraction, ...
        'windowStride',        opts.badchanstride, ...
        'averageReferenced',   opts.badchanavgref, ...
        'flatThresholdMicroV', opts.flatthreshold, ...
        'flatMaxBrokenTime',   opts.flatmaxbrokentime, ...
        'ramsaver',            opts.ramsaver);

    %%% Flat channels, on the data as recorded. Must come before the average reference:
    %%% a dead channel reads as one constant value here, but subtracting the common
    %%% average turns it into minus that average, which has real variance and is not
    %%% flat at all. Same window grid and duration rule as the correlation criterion.
    D = tic; fprintf('\nFlat channel detection ...\n')
    [flatmask, flatprop] = gedai.detectFlatChannels(EEG.data, EEG.srate, ...
        'windowseconds',  opts.badchanwindowseconds, ...
        'maxbrokentime',  opts.flatmaxbrokentime, ...
        'threshold',      opts.flatthreshold);
    KeepTime.FlatChannelDetection = toc(D);

    D = tic; fprintf('\nBad channel detection ...\n')

    %%% Strongly HP filter the data before RANSAC bad channel detection
    D = tic; fprintf('\nSlow activity removal ...\n')
    ICA_HiPassFilt_IIR = filterbank(EEG.srate, 'ICA_HiPassFilt_IIR');
    if opts.ramsaver
        % Overwrite each channel's row in place instead of filtering the whole
        % (transposed, double-cast) matrix at once - trades speed for lower peak memory.
        for iCh = 1:EEG.nbchan
            fprintf('Filtering channel %d/%d ...\n', iCh, EEG.nbchan);
            EEG.data(iCh,:) = single(filtfilt(ICA_HiPassFilt_IIR, double(EEG.data(iCh,:))'))';
        end
    else
        EEG.data = single(filtfilt(ICA_HiPassFilt_IIR, double(EEG.data)')');
    end
    KeepTime.HPFilter = toc(D);

    %%% Average-reference for detection only, then undo it (see the note in the module
    %%% help above): a single-electrode reference drives the ring of channels around it
    %%% towards a spuriously low reconstruction correlation.
    avgRef = [];
    if opts.badchanavgref
        fprintf('Average-referencing for detection ...\n')
        % avgRef = median([EEG.data; zeros(1, size(EEG.data,2))], 1);
        avgRef = sum(EEG.data, 1) / (size(EEG.data, 1) + 1);
        EEG.data = EEG.data - avgRef;
    end

    %%% The flat mask is folded in inside the cached call, so the mask on disk is the
    %%% one actually applied - bidsfun_gedai reads it back to index the leadfield.
    [removed_channels, corrs, znoise, flatprop] = smartcache( ...
        @() detectBadChannels(EEG, flatmask, flatprop, bcp), ...
        badchanFile, opts.refresh, ...
        {'', 'removed_channels', 'corr', 'znoise', 'flatprop'});

    %%% Recover channels that only failed the line-noise criterion: near-reference
    %%% vertex electrodes systematically read a spurious noise ratio (see the note in
    %%% the module help above), so a channel whose correlation performance is otherwise
    %%% fine is not actually bad, just close to the reference - it no longer counts
    %%% towards removed_channels, and is marked green (instead of red) in the noise
    %%% topoplot below so the override is visible rather than silent.
    savedchans = find(znoise' > opts.badchannoisethreshold & ...
        (sum(corrs < opts.badchancorrthreshold, 2) ./ size(corrs, 2)) < opts.badchanminbrokentime);
    removed_channels(savedchans) = false;

    if ~isempty(avgRef)
        EEG.data = EEG.data + avgRef;
        clear avgRef
    end

    EEG.etc.badchans = struct('mask', removed_channels, 'corr', corrs, ...
        'znoise', znoise, 'flatprop', flatprop, 'params', bcp, 'savedchans', savedchans);
    KeepTime.BadChannelDetection = toc(D);

    EEG.etc.filterparams.BadChannels = bcp;
    EEG.etc.filterparams.BadChannels.nRemoved = nnz(removed_channels);
    EEG.etc.filterparams.BadChannels.nRecovered = numel(savedchans);

    bc = EEG.etc.badchans;

    %%% Figures: where the bad channels are, and when they went bad
    gedai.plotBadChannels(bc.corr, bc.znoise, EEG.urchanlocs, ...
        fullfile(figDir, [fileID '_desc-' opts.desc '_BadChannelTopoplot.png']), ...
        bc.flatprop, bc.params, bc.savedchans);
    gedai.plotBadChannelTime(bc.corr, bc.mask, ...
        fullfile(figDir, [fileID '_desc-' opts.desc '_BadChannelTimecourse.png']), ...
        'windowseconds', bc.params.windowSeconds * bc.params.windowStride, ...
        'corrThreshold', bc.params.corrThreshold, ...
        'title', fileID);

    %%% channels.tsv: the standard's own way to carry channel status, and the file to
    %%% edit by hand when a call needs overruling. The .mat keeps what does not fit in a
    %%% table (the per-window correlation matrix).
    %%% status_description names why a channel was removed, so it lists only the criteria
    %%% that remove: correlation and flatness. Line noise is reported instead, in its own
    %%% column, because it no longer decides anything (see the note above).
    lowcorrprop = sum(bc.corr < bc.params.corrThreshold, 2) ./ size(bc.corr, 2);
    byCorr      = lowcorrprop(:) > bc.params.maxBrokenTime;
    byFlat      = bc.flatprop(:) > bc.params.flatMaxBrokenTime;
    reason      = repmat("good", numel(bc.mask), 1);
    reason(byCorr) = "low_correlation";
    reason(byFlat) = "flat";                        % most specific wins

    T = table(string({EEG.urchanlocs.labels})', ...
              repmat("eeg", numel(bc.mask), 1), ...
              repmat("microV", numel(bc.mask), 1), ...
              string(repmat("good", numel(bc.mask), 1)), ...
              reason, lowcorrprop(:), bc.znoise(:), bc.flatprop(:), ...
        'VariableNames', {'name', 'type', 'units', 'status', 'status_description', ...
                          'low_correlation_prop', 'znoise', 'flat_prop'});
    T.status(bc.mask(:)) = "bad";
    writetable(T, fullfile(outDir, [fileID '_desc-' opts.desc '_channels.tsv']), ...
        'FileType', 'text', 'Delimiter', '\t');

    %%% JSON sidecar: timings plus the parameters detection actually ran with
    prepParams = struct();
    if isfield(EEG.etc, 'filterparams'), prepParams = EEG.etc.filterparams; end
    prepParams.targetSampleRate = opts.targetsrate;
    prepParams.removeDC         = opts.removeDC;
    prepParams.zeroPatchSeconds = opts.zeropatchseconds;
    sidecarjson(KeepTime, fullfile(outDir, [fileID '_desc-' opts.desc '_badchans.json']), ...
        struct('BadChannelParameters', prepParams));

    fprintf('%d of %d channels flagged: %s\n', nnz(bc.mask), numel(bc.mask), ...
        mat2str(find(bc.mask(:))'))

    catch ME
        fprintf('[ERROR] %s: %s\n', fileID, ME.message);
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
    fid = fopen(fullfile(opts.savepath, 'failed_files_badchans.json'), 'w');
    fprintf(fid, '%s', jsonencode([failures{:}], 'PrettyPrint', true));
    fclose(fid);
end
end

% -------------------------------------------------------------------------
function [signal, removed_channels, corrs, znoise, flatprop] = detectBadChannels(EEG, flatmask, flatprop, bcp)
% Union of the clean_channels criteria and the flat-line criterion, as one cached unit.
% Kept together so the mask written to the cache is the mask actually applied to the
% data: bidsfun_gedai reads it back to pick the matching rows of the leadfield, and a
% cache holding only part of the criteria would silently desync from the saved file.
[signal, removed_channels, corrs, znoise] = clean_channels(EEG, ...
    bcp.corrThreshold, bcp.noiseThreshold, bcp.windowSeconds, bcp.maxBrokenTime, ...
    bcp.numSamples, bcp.subsetSizeFraction, bcp.windowStride, bcp.ramsaver);
removed_channels = removed_channels(:) | flatmask(:);
end
