function report = bidsfun_check_outputs(BIDS, opts)
% BIDSFUN_CHECK_OUTPUTS  Inventory the prep-pipeline outputs and flag what is missing.
%
%   report = bidsfun_check_outputs(BIDS, Name, Value, ...)
%
%   Walks the same EEG recordings that bidsfun_detect_badchans, bidsfun_hp_zap_cleanline,
%   bidsfun_evalfigs, bidsfun_gedai and ICA/run-pamica.py iterate over (same bids.query
%   on task/acq), and for each recording checks whether every file those stages are
%   supposed to write actually exists on disk. Prints a per-recording status matrix and
%   the list of files still to produce, and returns the whole thing as a struct for
%   scripting.
%
%   Nothing is (re)computed and nothing is written - this is a read-only status check.
%
%   Name-Value (defaults mirror PrepPipelineDROP.m)
%   ----------------------------------------------
%   derivfolder     derivatives sub-folder            (default 'prep-zc-ged')
%   badchandesc     desc of bidsfun_detect_badchans   (default 'hp')
%   filtdesc        desc of bidsfun_hp_zap_cleanline  (default 'hpzc')
%   geddesc         desc of bidsfun_gedai             (default 'hpzcged')
%   icadesc         OUT_DESC of run-pamica.py         (default 'pamica')
%   refresh         treat every expected output as still-to-produce, even when it is
%                   already on disk - i.e. preview what a pipeline run with refresh=true
%                   would (re)generate. Present files are then marked STALE rather than
%                   ok, and all of them are listed under "still to produce". (default false)
%   tasklabel       BIDS task label(s)                (default {'Sleep','sleep'})
%   acqlabel        BIDS acq label                    (default '')
%   subjectfilter   cell of subject IDs; {} = all
%   sessionfilter   cell of session IDs; {} = all
%   savepath        derivatives root  (default <BIDS.pth>/derivatives/<derivfolder>)
%   figpath         figures root      (default <savepath>/figures)
%   checkfigures    also check the evalfigs sentinel PNGs (default true)
%   checkica        also check the run-pamica.py / ICLabel outputs (default true)
%   csvout          path to also write the status table as CSV; '' = skip (default '')

arguments
    BIDS
    opts.derivfolder    char = 'prep-zc-ged'
    opts.badchandesc    char = 'hp'
    opts.filtdesc       char = 'hpzc'
    opts.geddesc        char = 'hpzcged'
    opts.icadesc        char = 'pamica'
    opts.refresh (1,1) logical = false
    opts.tasklabel           = {'Sleep', 'sleep'}
    opts.acqlabel       char = ''
    opts.subjectfilter  cell = {}
    opts.sessionfilter  cell = {}
    opts.savepath       char = ''
    opts.figpath        char = ''
    opts.checkfigures (1,1) logical = true
    opts.checkica     (1,1) logical = true
    opts.csvout         char = ''
end

if isempty(opts.savepath), opts.savepath = fullfile(BIDS.pth, 'derivatives', opts.derivfolder); end
if isempty(opts.figpath),  opts.figpath  = fullfile(opts.savepath, 'figures'); end

fprintf('\n=== bidsfun_check_outputs ===\n');
fprintf('derivatives : %s\n', opts.savepath);
fprintf('desc labels : badchan=%s  filt=%s  ged=%s  ica=%s\n', ...
    opts.badchandesc, opts.filtdesc, opts.geddesc, opts.icadesc);
if opts.refresh
    fprintf('refresh     : ON - every expected output is reported as still-to-produce\n');
end
fprintf('\n');

%%% Same recording query the pipeline stages use
q = {'data', 'extension', '.vhdr', 'task', opts.tasklabel};
if ~isempty(opts.acqlabel), q = [q, {'acq', opts.acqlabel}]; end
filesEEG = bids.query(BIDS, q{:});
if isempty(filesEEG)
    error('bidsfun_check_outputs:noFiles', 'No matching EEG files found in BIDS layout.');
end

%%% Column definitions: {label, stage, group, showInTable}
%%%   group: 'data' always checked, 'fig' gated by checkfigures, 'ica' gated by checkica
%%%   showInTable: sidecars/tsv are still checked and still appear in the completeness
%%%   summary and the "still to produce" list, but the per-recording matrix shows only
%%%   one representative file per stage to stay readable.
cols = { ...
    'badchan.mat'    1 'data' true
    'badchan.tsv'    1 'data' false
    'badchan.json'   1 'data' false
    'zc.set'         2 'data' true
    'zc.json'        2 'data' false
    'eval_zc.png'    2 'fig'  true
    'ged.set'        3 'data' true
    'eval_ged.png'   3 'fig'  true
    'ica.mat'        4 'ica'  true
    'ica.json'       4 'ica'  false
    'iclabels.tsv'   4 'ica'  true
    'iclabels.json'  4 'ica'  false };
colLabels = cols(:,1);
colGroup  = cols(:,3);
colShow   = cell2mat(cols(:,4))';
nCol      = size(cols,1);

active = true(1, nCol);
active(strcmp(colGroup, 'fig')) = opts.checkfigures;
active(strcmp(colGroup, 'ica')) = opts.checkica;

rows     = {};   % one struct per recording
todoList = {};   % flat list of files still to produce (respects refresh)
missList = {};   % flat list of files genuinely absent (ignores refresh)

for ifile = 1:numel(filesEEG)
    eegFile = filesEEG{ifile};
    p       = bids.internal.parse_filename(eegFile);
    fileID  = strjoin(cellfun(@(k) [k '-' p.entities.(k)], fieldnames(p.entities), 'uni', 0), '_');
    subDir  = fullfile(['sub-' p.entities.sub], ['ses-' p.entities.ses]);

    if ~isempty(opts.subjectfilter) && ~contains(fileID, opts.subjectfilter), continue, end
    if ~isempty(opts.sessionfilter) && ~contains(fileID, opts.sessionfilter), continue, end

    outDir = fullfile(opts.savepath, subDir);

    %%% run-pamica.py swaps desc-<geddesc> -> desc-<icadesc> and _eeg -> _ica
    icaBase = fullfile(outDir, sprintf('%s_desc-%s_ica', fileID, opts.icadesc));

    exp = containers.Map('KeyType', 'char', 'ValueType', 'char');
    exp('badchan.mat')   = fullfile(outDir, sprintf('%s_desc-%s_badchans.mat',  fileID, opts.badchandesc));
    exp('badchan.tsv')   = fullfile(outDir, sprintf('%s_desc-%s_channels.tsv',  fileID, opts.badchandesc));
    exp('badchan.json')  = fullfile(outDir, sprintf('%s_desc-%s_badchans.json', fileID, opts.badchandesc));
    exp('zc.set')        = fullfile(outDir, sprintf('%s_desc-%s_eeg.set',       fileID, opts.filtdesc));
    exp('zc.json')       = fullfile(outDir, sprintf('%s_desc-%s_eeg.json',      fileID, opts.filtdesc));
    exp('ged.set')       = fullfile(outDir, sprintf('%s_desc-%s_eeg.set',       fileID, opts.geddesc));
    exp('eval_zc.png')   = fullfile(opts.figpath, ['desc-' opts.filtdesc], subDir, [fileID '_psd_per_stage.png']);
    exp('eval_ged.png')  = fullfile(opts.figpath, ['desc-' opts.geddesc],  subDir, [fileID '_psd_per_stage.png']);
    exp('ica.mat')       = [icaBase '.mat'];
    exp('ica.json')      = [icaBase '.json'];
    exp('iclabels.tsv')  = fullfile(outDir, sprintf('%s_desc-%s_iclabels.tsv',  fileID, opts.icadesc));
    exp('iclabels.json') = fullfile(outDir, sprintf('%s_desc-%s_iclabels.json', fileID, opts.icadesc));

    r = struct('fileID', fileID, 'sub', ['sub-' p.entities.sub], 'ses', ['ses-' p.entities.ses]);
    status = nan(1, nCol);   % 1 present, 0 missing, 2 present-but-stale (refresh), NaN not checked
    for c = 1:nCol
        if ~active(c), continue, end
        f  = exp(colLabels{c});
        ok = isfile(f);
        if ~ok
            status(c) = 0;
            missList{end+1,1} = f; %#ok<AGROW>
            todoList{end+1,1} = f; %#ok<AGROW>
        elseif opts.refresh
            status(c) = 2;
            todoList{end+1,1} = f; %#ok<AGROW>
        else
            status(c) = 1;
        end
    end
    r.status = status;
    rows{end+1,1} = r; %#ok<AGROW>
end

if isempty(rows)
    fprintf('No recordings matched the subject/session filter.\n');
    report = struct('rows', {rows}, 'todo', {todoList}, 'missing', {missList}, 'colLabels', {colLabels});
    return
end

%%% ---- Per-recording status matrix ----
w   = max(cellfun(@(x) numel(x.fileID), rows)) + 2;
cw  = max(cellfun(@numel, colLabels)) + 2;
fmt = ['%-' num2str(cw) 's'];
showIdx = find(colShow & active);
fprintf(['%-' num2str(w) 's'], '');
for c = showIdx, fprintf(fmt, colLabels{c}); end
fprintf('\n');
for i = 1:numel(rows)
    r = rows{i};
    fprintf(['%-' num2str(w) 's'], r.fileID);
    for c = showIdx
        switch double(r.status(c))
            case 1,    mark = 'ok';
            case 2,    mark = 'STALE';
            case 0,    mark = 'MISSING';
            otherwise, mark = '-';
        end
        fprintf(fmt, mark);
    end
    fprintf('\n');
end

%%% ---- Completeness by output ----
fprintf('\n--- Completeness by output ---\n');
for c = 1:nCol
    if ~active(c), continue, end
    s = cellfun(@(r) r.status(c), rows);
    fprintf('  %-14s %3d / %3d present\n', colLabels{c}, sum(s == 1), numel(s));
end

%%% ---- Files still to produce ----
label = 'Missing files';
if opts.refresh, label = 'Files a refresh run would (re)produce'; end
fprintf('\n--- %s (%d) ---\n', label, numel(todoList));
for i = 1:numel(todoList)
    fprintf('  %s\n', todoList{i});
end
if isempty(todoList)
    fprintf('  none - every expected output is present.\n');
end
if opts.refresh && ~isempty(missList)
    fprintf('\n  (of which %d are genuinely absent, the rest would be overwritten)\n', numel(missList));
end

%%% ---- Table / optional CSV ----
%%% Same columns as the printed matrix - one representative file per stage.
statusMat = cell2mat(cellfun(@(r) r.status, rows, 'uni', 0));
T = array2table(statusMat(:, showIdx), 'VariableNames', matlab.lang.makeValidName(colLabels(showIdx)));
T = addvars(T, string(cellfun(@(r) r.fileID, rows, 'uni', 0)), 'Before', 1, 'NewVariableNames', 'fileID');
if ~isempty(opts.csvout)
    writetable(T, opts.csvout);
    fprintf('\nStatus table written to %s\n', opts.csvout);
end

report = struct('table', T, 'rows', {rows}, 'todo', {todoList}, 'missing', {missList}, ...
    'colLabels', {colLabels}, 'files', {filesEEG});
end
