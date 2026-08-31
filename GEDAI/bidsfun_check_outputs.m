function report = bidsfun_check_outputs(BIDS, opts)
% BIDSFUN_CHECK_OUTPUTS  Inventory the outputs of the prep pipeline and flag what is missing.
%
%   report = bidsfun_check_outputs(BIDS, Name, Value, ...)
%
%   Walks the same EEG recordings that bidsfun_detect_badchans, bidsfun_hp_zap_cleanline,
%   bidsfun_evalfigs and bidsfun_gedai iterate over (same bids.query on task/acq), and for
%   each recording checks whether every file those four stages are supposed to write
%   actually exists on disk. Prints a per-recording status matrix and a list of the
%   missing files, and returns the whole thing as a struct/table for scripting.
%
%   Nothing is (re)computed and nothing is written - this is a read-only status check.
%
%   Name-Value (defaults mirror PrepPipelineDROP.m)
%   ----------------------------------------------
%   derivfolder     derivatives sub-folder            (default 'prep-zc-ged')
%   badchandesc     desc of bidsfun_detect_badchans   (default 'hp')
%   filtdesc        desc of bidsfun_hp_zap_cleanline  (default 'hpzc')
%   geddesc         desc of bidsfun_gedai             (default 'hpzcged')
%   tasklabel       BIDS task label(s)                (default {'Sleep','sleep'})
%   acqlabel        BIDS acq label                    (default '')
%   subjectfilter   cell of subject IDs; {} = all
%   sessionfilter   cell of session IDs; {} = all
%   savepath        derivatives root  (default <BIDS.pth>/derivatives/<derivfolder>)
%   figpath         figures root      (default <savepath>/figures)
%   checkfigures    also check the evalfigs sentinel PNGs (default true)
%   csvout          path to also write the status table as CSV; '' = skip (default '')

arguments
    BIDS
    opts.derivfolder    char = 'prep-zc-ged'
    opts.badchandesc    char = 'hp'
    opts.filtdesc       char = 'hpzc'
    opts.geddesc        char = 'hpzcged'
    opts.tasklabel           = {'Sleep', 'sleep'}
    opts.acqlabel       char = ''
    opts.subjectfilter  cell = {}
    opts.sessionfilter  cell = {}
    opts.savepath       char = ''
    opts.figpath        char = ''
    opts.checkfigures (1,1) logical = true
    opts.csvout         char = ''
end

if isempty(opts.savepath), opts.savepath = fullfile(BIDS.pth, 'derivatives', opts.derivfolder); end
if isempty(opts.figpath),  opts.figpath  = fullfile(opts.savepath, 'figures'); end

fprintf('\n=== bidsfun_check_outputs ===\n');
fprintf('derivatives : %s\n', opts.savepath);
fprintf('desc labels : badchan=%s  filt=%s  ged=%s\n\n', ...
    opts.badchandesc, opts.filtdesc, opts.geddesc);

%%% Same recording query the four stages use
if isempty(opts.acqlabel)
    filesEEG = bids.query(BIDS, 'data', 'extension', '.vhdr', 'task', opts.tasklabel);
else
    filesEEG = bids.query(BIDS, 'data', 'extension', '.vhdr', 'task', opts.tasklabel, 'acq', opts.acqlabel);
end
if isempty(filesEEG)
    error('bidsfun_check_outputs:noFiles', 'No matching EEG files found in BIDS layout.');
end

%%% Column definitions: {label, stage, isFigure}
cols = { ...
    'badchan.mat'   1 false
    'badchan.tsv'   1 false
    'badchan.json'  1 false
    'filt.set'      2 false
    'filt.json'     2 false
    'evalfilt.png'  2 true
    'ged.set'       3 false
    'ged.json'      3 false
    'evalged.png'   3 true };
colLabels = cols(:,1);
nCol      = size(cols,1);

rows        = {};   % one struct per recording
missingList = {};   % flat list of missing file paths

for ifile = 1:numel(filesEEG)
    eegFile = filesEEG{ifile};
    p       = bids.internal.parse_filename(eegFile);
    fileID  = strjoin(cellfun(@(k) [k '-' p.entities.(k)], fieldnames(p.entities), 'uni', 0), '_');
    subDir  = fullfile(['sub-' p.entities.sub], ['ses-' p.entities.ses]);

    if ~isempty(opts.subjectfilter) && ~contains(fileID, opts.subjectfilter), continue, end
    if ~isempty(opts.sessionfilter) && ~contains(fileID, opts.sessionfilter), continue, end

    outDir = fullfile(opts.savepath, subDir);

    %%% Expected files, keyed to colLabels above
    exp = containers.Map('KeyType', 'char', 'ValueType', 'char');
    exp('badchan.mat')  = fullfile(outDir, sprintf('%s_desc-%s_badchans.mat',  fileID, opts.badchandesc));
    exp('badchan.tsv')  = fullfile(outDir, sprintf('%s_desc-%s_channels.tsv',  fileID, opts.badchandesc));
    exp('badchan.json') = fullfile(outDir, sprintf('%s_desc-%s_badchans.json', fileID, opts.badchandesc));
    exp('filt.set')     = fullfile(outDir, sprintf('%s_desc-%s_eeg.set',       fileID, opts.filtdesc));
    exp('filt.json')    = fullfile(outDir, sprintf('%s_desc-%s_eeg.json',      fileID, opts.filtdesc));
    exp('ged.set')      = fullfile(outDir, sprintf('%s_desc-%s_eeg.set',       fileID, opts.geddesc));
    exp('ged.json')     = fullfile(outDir, sprintf('%s_desc-%s_eeg.json',      fileID, opts.geddesc));
    exp('evalfilt.png') = fullfile(opts.figpath, ['desc-' opts.filtdesc], subDir, [fileID '_psd_per_stage.png']);
    exp('evalged.png')  = fullfile(opts.figpath, ['desc-' opts.geddesc],  subDir, [fileID '_psd_per_stage.png']);

    r = struct('fileID', fileID, 'sub', ['sub-' p.entities.sub], 'ses', ['ses-' p.entities.ses]);
    status = zeros(1, nCol);   % 1 = present, 0 = missing, NaN = not checked
    for c = 1:nCol
        key = colLabels{c};
        if cols{c,3} && ~opts.checkfigures
            status(c) = NaN;         % figure check disabled -> not counted
            continue
        end
        ok = isfile(exp(key));
        status(c) = ok;
        if ~ok
            missingList{end+1,1} = exp(key); %#ok<AGROW>
        end
    end
    r.status  = status;
    r.missing = colLabels(status == 0);
    rows{end+1,1} = r; %#ok<AGROW>
end

if isempty(rows)
    fprintf('No recordings matched the subject/session filter.\n');
    report = struct('rows', {rows}, 'missing', {missingList}, 'colLabels', {colLabels});
    return
end

%%% ---- Per-recording status matrix ----
w = max(cellfun(@(x) numel(x.fileID), rows)) + 2;
hdr = [repmat(' ',1,w) sprintf('%-14s', colLabels{:})];
fprintf('%s\n', hdr);
for i = 1:numel(rows)
    r = rows{i};
    line = sprintf(['%-' num2str(w) 's'], r.fileID);
    for c = 1:nCol
        switch double(r.status(c))
            case 1, mark = 'ok';
            case 0, mark = 'MISSING';
            otherwise, mark = '-';
        end
        line = [line sprintf('%-14s', mark)]; %#ok<AGROW>
    end
    fprintf('%s\n', line);
end

%%% ---- Per-stage completeness ----
fprintf('\n--- Completeness by output ---\n');
for c = 1:nCol
    s = cellfun(@(r) r.status(c), rows);
    if all(isnan(s)), continue, end
    fprintf('  %-14s %3d / %3d present\n', colLabels{c}, sum(s == 1), sum(~isnan(s)));
end

%%% ---- Missing files ----
fprintf('\n--- Missing files (%d) ---\n', numel(missingList));
for i = 1:numel(missingList)
    fprintf('  %s\n', missingList{i});
end
if isempty(missingList)
    fprintf('  none - every expected output is present.\n');
end

%%% ---- Optional CSV ----
statusMat = cell2mat(cellfun(@(r) r.status, rows, 'uni', 0));
T = array2table(statusMat, 'VariableNames', matlab.lang.makeValidName(colLabels));
T = addvars(T, string(cellfun(@(r) r.fileID, rows, 'uni', 0)), 'Before', 1, 'NewVariableNames', 'fileID');
if ~isempty(opts.csvout)
    writetable(T, opts.csvout);
    fprintf('\nStatus table written to %s\n', opts.csvout);
end

report = struct('table', T, 'rows', {rows}, 'missing', {missingList}, ...
    'colLabels', {colLabels}, 'files', {filesEEG});
end
