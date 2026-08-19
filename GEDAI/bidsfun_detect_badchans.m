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
%   badchanavgref   average-reference for detection              (default true)
%   badchanstride   evaluate every Nth window in the correlation criterion (default 2)
%   flatthreshold   peak-to-peak uV below which a window is flat (default 0.5)
%   noteegchannels  channel indices to drop                      (default 257:300)
%   sfppath         path passed to the SFP resolver              (default BIDS root)
%   subjectfilter   cell array of subject ID strings; {} = all subjects
%   sessionfilter   cell array of session ID strings; {} = all sessions
%
%   See run.run_filter for what the criteria are and why they sit where they do.

arguments
    BIDS

    %--- Paths ---
    opts.derivfolder      char    = 'prep-ged'
    opts.savepath         char    = ''
    opts.figpath          char    = ''
    opts.refresh (1,1)    logical = false
    opts.desc             char    = 'badchan'

    %--- EEG ---
    opts.tasklabel                       = {'Sleep', 'sleep'}
    opts.acqlabel   char                 = ''
    opts.noteegchannels   (1,:) double   = 257:300
    opts.targetsrate      (1,1) double   = 0
    opts.removeDC         (1,1) logical  = true
    opts.zeropatchseconds (1,1) double   = 5

    %--- Bad channels ---
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

    %%% Detect. run.run_filter stops after detection under badchannelsonly, so the
    %%% preparation is identical to what bidsfun_hp_zap_cleanline will apply later.
    [EEG, KeepTime] = run.run_filter(EEG, ...
        'noteegchannels',     opts.noteegchannels, ...
        'targetsrate',        opts.targetsrate, ...
        'removeDC',           opts.removeDC, ...
        'zeropatchseconds',   opts.zeropatchseconds, ...
        'restorezeropatches', false, ...
        'KeepTime',           KeepTime, ...
        'badchannels',        true, ...
        'badchannelsonly',    true, ...
        'badchanavgref',      opts.badchanavgref, ...
        'badchanstride',      opts.badchanstride, ...
        'flatthreshold',      opts.flatthreshold, ...
        'badchanfile',        badchanFile, ...
        'refresh',            opts.refresh);

    bc = EEG.etc.badchans;

    %%% Figures: where the bad channels are, and when they went bad
    gedai.plotBadChannels(bc.corr, bc.znoise, EEG.urchanlocs, ...
        fullfile(figDir, [fileID '_desc-' opts.desc '_BadChannelTopoplot.png']), ...
        bc.flatprop, bc.params);
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
    %%% column, because it no longer decides anything (see the note in run.run_filter).
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
