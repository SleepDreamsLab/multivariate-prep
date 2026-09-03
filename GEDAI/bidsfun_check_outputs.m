function report = bidsfun_check_outputs(BIDS, opts)
% BIDSFUN_CHECK_OUTPUTS  Inventory the prep-pipeline inputs and outputs, and flag what is missing.
%
%   report = bidsfun_check_outputs(BIDS, Name, Value, ...)
%
%   Walks the same EEG recordings that bidsfun_detect_badchans, bidsfun_hp_zap_cleanline,
%   bidsfun_evalfigs, bidsfun_gedai and ICA/run-pamica.py iterate over (same bids.query
%   on task/acq), and per recording checks
%     - the external INPUTS bidsfun_gedai needs: sleep scoring, SFP montage, leadfield
%     - the OUTPUTS every stage is supposed to write
%   Prints a per-recording status matrix, a completeness summary and the list of files
%   still to produce, draws a heatmap, and returns the whole thing as a struct.
%
%   The two are worth separating: a missing output just means the stage has not run yet,
%   while a missing input means it can never run for that recording until someone
%   produces the file - so they are listed apart, and the heatmap rules them off.
%
%   Inputs are resolved with the pipeline's own resolvers (gedai.collectScoringFiles +
%   gedai.matchScoringFile, gedai.matchSfpFile, and gedai.loadrefcov's path rule) rather
%   than a reimplementation, so this check cannot drift from what bidsfun_gedai does.
%
%   Nothing is (re)computed and nothing is written except the heatmap/CSV.
%
%   Name-Value (defaults mirror PrepPipelineDROP.m / bidsfun_gedai)
%   --------------------------------------------------------------
%   derivfolder     derivatives sub-folder            (default 'prep-zc-ged')
%   badchandesc     desc of bidsfun_detect_badchans   (default 'hp')
%   filtdesc        desc of bidsfun_hp_zap_cleanline  (default 'hpzc')
%   geddesc         desc of bidsfun_gedai             (default 'hpzcged')
%   icadesc         OUT_DESC of run-pamica.py         (default 'pamica')
%   refresh         treat every expected OUTPUT as still-to-produce, even when it is
%                   already on disk - i.e. preview what a pipeline run with refresh=true
%                   would regenerate. Inputs are unaffected: nothing here produces them.
%                   (default false)
%   tasklabel       BIDS task label(s)                (default {'Sleep','sleep'})
%   acqlabel        BIDS acq label                    (default '')
%   subjectfilter   cell of subject IDs; {} = all
%   sessionfilter   cell of session IDs; {} = all
%   savepath        derivatives root  (default <BIDS.pth>/derivatives/<derivfolder>)
%   figpath         figures root      (default <savepath>/figures)
%   checkinputs     check scoring / SFP / leadfield   (default true)
%   scoringpath     scoring root, as passed to bidsfun_gedai. '' searches
%                   <BIDS.pth>/<sub>/<ses> per recording, exactly as that function does.
%                   (default <BIDS.pth>/derivatives/scoring/scores/Manual_Checked)
%   sfppath         path handed to the SFP resolver   (default BIDS.pth)
%   leadfieldpath   leadfield root; a folder is read as
%                   <root>/<sub>/<ses>/headmodel_surf_openmeeg.mat, anything else as the
%                   file itself
%                   (default <BIDS.pth>/../Data_Analysis/Brainstorm_db/Leadfield_PM/data)
%   checkfigures    also check the evalfigs sentinel PNGs (default true)
%   checkica        also check the run-pamica.py / ICLabel outputs (default true)
%   csvout          path to also write the status table as CSV; '' = skip (default '')
%   plot            draw the status heatmap             (default true)
%   plotfile        where to save it ('' = <figpath>/prep_status_heatmap.png)

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

    %--- External inputs (same defaults as bidsfun_gedai) ---
    opts.checkinputs (1,1) logical = true
    opts.scoringpath    char = fullfile(BIDS.pth, 'derivatives', 'scoring', 'scores', 'Manual_Checked')
    opts.sfppath        char = BIDS.pth
    opts.leadfieldpath  char = fullfile(BIDS.pth, '..', 'Data_Analysis', 'Brainstorm_db', 'Leadfield_PM', 'data')

    opts.checkfigures (1,1) logical = true
    opts.checkica     (1,1) logical = true
    opts.csvout         char = ''
    opts.plot   (1,1) logical = true
    opts.plotfile       char = ''
end

if isempty(opts.savepath), opts.savepath = fullfile(BIDS.pth, 'derivatives', opts.derivfolder); end
if isempty(opts.figpath),  opts.figpath  = fullfile(opts.savepath, 'figures'); end

fprintf('\n=== bidsfun_check_outputs ===\n');
fprintf('derivatives : %s\n', opts.savepath);
fprintf('desc labels : badchan=%s  filt=%s  ged=%s  ica=%s\n', ...
    opts.badchandesc, opts.filtdesc, opts.geddesc, opts.icadesc);
if opts.checkinputs
    fprintf('scoring     : %s\n', ternary(isempty(opts.scoringpath), '<per recording, beside the raw data>', opts.scoringpath));
    fprintf('sfp         : %s\n', opts.sfppath);
    fprintf('leadfields  : %s\n', opts.leadfieldpath);
end
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

%%% Scoring files, collected once for the whole tree - the same single recursive scan
%%% bidsfun_gedai does when scoringpath is set. With scoringpath empty it collects per
%%% recording instead (below), again mirroring that function.
scoringfiles = {};
if opts.checkinputs && ~isempty(opts.scoringpath)
    if ~isfolder(opts.scoringpath)
        warning('bidsfun_check_outputs:noScoringDir', ...
            'Scoring directory does not exist: %s - every recording will report scoring as missing.', ...
            opts.scoringpath);
    else
        fprintf('Scanning scoring directory ...\n');
        scoringfiles = gedai.collectScoringFiles(opts.scoringpath);
        fprintf('  %d scoring file(s) found\n\n', numel(scoringfiles));
    end
end

%%% gedai.matchSfpFile warns once per ambiguous session; over a few hundred recordings
%%% that buries the actual report, and the ambiguity is a naming problem for the montage
%%% folder rather than something this check can act on.
wstate = warning('off', 'gedai:matchSfpFile:ambiguous');
cleanupWarn = onCleanup(@() warning(wstate));

%%% Column definitions: {label, stage, group, showInTable}
%%%   group: 'input' gated by checkinputs, 'data' always checked,
%%%          'fig' gated by checkfigures, 'ica' gated by checkica
%%%   showInTable: sidecars/tsv are still checked and still appear in the completeness
%%%   summary and the "still to produce" list, but the per-recording matrix shows only
%%%   one representative file per stage to stay readable.
cols = { ...
    'scoring'        0 'input' true
    'sfp'            0 'input' true
    'leadfield'      0 'input' true
    'badchan.mat'    1 'data'  true
    'badchan.tsv'    1 'data'  false
    'badchan.json'   1 'data'  false
    'zc.set'         2 'data'  true
    'zc.json'        2 'data'  false
    'eval_zc.png'    2 'fig'   true
    'ged.set'        3 'data'  true
    'eval_ged.png'   3 'fig'   true
    'ica.mat'        4 'ica'   true
    'ica.json'       4 'ica'   false
    'iclabels.tsv'   4 'ica'   true
    'iclabels.json'  4 'ica'   false };
colLabels = cols(:,1);
colGroup  = cols(:,3);
colShow   = cell2mat(cols(:,4))';
nCol      = size(cols,1);

isInput = strcmp(colGroup, 'input')';
active  = true(1, nCol);
active(isInput)                 = opts.checkinputs;
active(strcmp(colGroup, 'fig')) = opts.checkfigures;
active(strcmp(colGroup, 'ica')) = opts.checkica;

rows      = {};   % one struct per recording
todoList  = {};   % outputs still to produce (respects refresh)
missList  = {};   % outputs genuinely absent (ignores refresh)
blockList = {};   % missing INPUTS - nothing in the pipeline produces these

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

    %%% Expected outputs: the path is known up front, so presence is one isfile each.
    expect = containers.Map('KeyType', 'char', 'ValueType', 'char');
    found  = containers.Map('KeyType', 'char', 'ValueType', 'logical');
    expect('badchan.mat')   = fullfile(outDir, sprintf('%s_desc-%s_badchans.mat',  fileID, opts.badchandesc));
    expect('badchan.tsv')   = fullfile(outDir, sprintf('%s_desc-%s_channels.tsv',  fileID, opts.badchandesc));
    expect('badchan.json')  = fullfile(outDir, sprintf('%s_desc-%s_badchans.json', fileID, opts.badchandesc));
    expect('zc.set')        = fullfile(outDir, sprintf('%s_desc-%s_eeg.set',       fileID, opts.filtdesc));
    expect('zc.json')       = fullfile(outDir, sprintf('%s_desc-%s_eeg.json',      fileID, opts.filtdesc));
    expect('ged.set')       = fullfile(outDir, sprintf('%s_desc-%s_eeg.set',       fileID, opts.geddesc));
    expect('eval_zc.png')   = fullfile(opts.figpath, ['desc-' opts.filtdesc], subDir, [fileID '_psd_per_stage.png']);
    expect('eval_ged.png')  = fullfile(opts.figpath, ['desc-' opts.geddesc],  subDir, [fileID '_psd_per_stage.png']);
    expect('ica.mat')       = [icaBase '.mat'];
    expect('ica.json')      = [icaBase '.json'];
    expect('iclabels.tsv')  = fullfile(outDir, sprintf('%s_desc-%s_iclabels.tsv',  fileID, opts.icadesc));
    expect('iclabels.json') = fullfile(outDir, sprintf('%s_desc-%s_iclabels.json', fileID, opts.icadesc));
    for c = find(~isInput)
        found(colLabels{c}) = isfile(expect(colLabels{c}));
    end

    %%% External inputs: resolved by search, so the path is only known once found. When
    %%% nothing matches, what gets recorded is where it looked - that is the actionable
    %%% part of a missing input.
    if opts.checkinputs
        [expect('scoring'),   found('scoring')]   = resolveScoring(p, subDir, BIDS, opts, scoringfiles);
        [expect('sfp'),       found('sfp')]       = resolveSfp(p, opts);
        [expect('leadfield'), found('leadfield')] = resolveLeadfield(p, opts);
    end

    r = struct('fileID', fileID, 'sub', ['sub-' p.entities.sub], 'ses', ['ses-' p.entities.ses]);
    status = nan(1, nCol);   % 1 present, 0 missing, 2 present-but-stale (refresh), NaN not checked
    for c = 1:nCol
        if ~active(c), continue, end
        key = colLabels{c};
        f   = expect(key);
        if ~found(key)
            status(c) = 0;
            if isInput(c)
                blockList{end+1,1} = sprintf('%-42s %s : %s', fileID, key, f); %#ok<AGROW>
            else
                missList{end+1,1} = f; %#ok<AGROW>
                todoList{end+1,1} = f; %#ok<AGROW>
            end
        elseif opts.refresh && ~isInput(c)
            %%% Inputs are never regenerated, so refresh does not apply to them.
            status(c) = 2;
            todoList{end+1,1} = f; %#ok<AGROW>
        else
            status(c) = 1;
        end
    end
    r.status = status;
    if opts.checkinputs
        r.inputs = struct('scoring', expect('scoring'), 'sfp', expect('sfp'), ...
            'leadfield', expect('leadfield'));
    end
    rows{end+1,1} = r; %#ok<AGROW>
end

if isempty(rows)
    fprintf('No recordings matched the subject/session filter.\n');
    report = struct('rows', {rows}, 'todo', {todoList}, 'missing', {missList}, ...
        'blockers', {blockList}, 'colLabels', {colLabels});
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

%%% ---- Completeness ----
fprintf('\n--- Completeness by file ---\n');
for c = 1:nCol
    if ~active(c), continue, end
    s = cellfun(@(r) r.status(c), rows);
    fprintf('  %-14s %3d / %3d present%s\n', colLabels{c}, sum(s ~= 0), numel(s), ...
        ternary(isInput(c), '   (input)', ''));
end

%%% ---- Missing inputs ----
if opts.checkinputs
    fprintf('\n--- Missing INPUTS (%d) - these block the stages that need them ---\n', numel(blockList));
    for i = 1:numel(blockList)
        fprintf('  %s\n', blockList{i});
    end
    if isempty(blockList)
        fprintf('  none - scoring, montage and leadfield resolve for every recording.\n');
    end
end

%%% ---- Outputs still to produce ----
label = 'Missing OUTPUTS';
if opts.refresh, label = 'OUTPUTS a refresh run would (re)produce'; end
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

%%% ---- Heatmap ----
fig = [];
if opts.plot
    if isempty(opts.plotfile)
        opts.plotfile = fullfile(opts.figpath, 'prep_status_heatmap.png');
    end
    fig = plotStatusHeatmap(statusMat(:, showIdx), colLabels(showIdx), ...
        cellfun(@(r) r.fileID, rows, 'uni', 0), opts.plotfile, ...
        nnz(isInput(showIdx)));
end

report = struct('table', T, 'rows', {rows}, 'todo', {todoList}, 'missing', {missList}, ...
    'blockers', {blockList}, 'colLabels', {colLabels}, 'files', {filesEEG}, 'figure', fig);
end

% -------------------------------------------------------------------------
function [pathOrWhere, ok] = resolveScoring(p, subDir, BIDS, opts, scoringfiles)
% Sleep scoring, matched exactly as bidsfun_gedai does: collect the candidates, then
% match on every entity in the EEG filename except acq. An empty scoringpath means the
% scores sit beside the raw recording, so the collection happens per recording there.
where = opts.scoringpath;
if isempty(opts.scoringpath)
    where = fullfile(BIDS.pth, subDir);
    scoringfiles = gedai.collectScoringFiles(where);
end
f = gedai.matchScoringFile(p.entities, scoringfiles);
ok = ~isempty(f);
pathOrWhere = ternary(ok, f, sprintf('no scoring matched under %s', where));
end

% -------------------------------------------------------------------------
function [pathOrWhere, ok] = resolveSfp(p, opts)
% SFP montage, via the same study-dispatching resolver gedai.assignChanlocs calls. That
% resolver errors outright on a subject name it cannot map to a study - caught here, since
% one unknown study should report as a missing montage rather than abort the whole scan.
try
    f  = gedai.matchSfpFile(opts.sfppath, p.entities.sub, p.entities.ses);
    ok = ~isempty(f);
    pathOrWhere = ternary(ok, f, sprintf('no SFP matched for sub-%s ses-%s under %s', ...
        p.entities.sub, p.entities.ses, opts.sfppath));
catch ME
    ok = false;
    pathOrWhere = sprintf('SFP resolver failed: %s', ME.message);
end
end

% -------------------------------------------------------------------------
function [pathOrWhere, ok] = resolveLeadfield(p, opts)
% Brainstorm headmodel, addressed the way gedai.loadrefcov addresses it: a folder is read
% as <root>/sub-X/ses-Y/headmodel_surf_openmeeg.mat, anything else as the file itself.
if isfolder(opts.leadfieldpath)
    f = fullfile(opts.leadfieldpath, ['sub-' p.entities.sub], ['ses-' p.entities.ses], ...
        'headmodel_surf_openmeeg.mat');
else
    f = opts.leadfieldpath;
end
ok = isfile(f);
pathOrWhere = f;
end

% -------------------------------------------------------------------------
function fig = plotStatusHeatmap(S, colLabels, fileIDs, savefile, nInputCols)
% imagesc grid: rows = recordings, columns = one file per pipeline stage, with the
% external inputs ruled off from the outputs the pipeline itself writes.
%   0 missing (red) | 1 present (green) | NaN not checked (grey)
nRow = size(S, 1);
M = S;
M(M == 2)   = 1;                 % "stale" (refresh) shown the same as present
M(isnan(M)) = -1;                % fold "not checked" into its own colour bin

%%% Rows are grouped by subject, so a participant label heads a BLOCK of recordings, not a
%%% single row - both the tick placement and the figure height follow from that.
subj     = regexprep(fileIDs(:), '_.*$', '');
firstIdx = find([true; ~strcmp(subj(2:end), subj(1:end-1))]);

fig = figure('Color', 'w', 'Name', 'prep status', ...
    'Position', [100 100 max(560, 90*numel(colLabels)+260) ...
                 min(1200, 260 + max(14*nRow, 16*numel(firstIdx)))]);
ax = axes(fig);
imagesc(ax, M);
cmap = [0.75 0.75 0.75;   % -1 not checked
        0.85 0.20 0.20;   %  0 missing
        0.20 0.65 0.30];  %  1 present
colormap(ax, cmap);
set(ax, 'CLim', [-1.5 1.5]);

ax.XTick = 1:numel(colLabels);
ax.XTickLabel = colLabels;
ax.XTickLabelRotation = 30;
ax.TickLabelInterpreter = 'none';
ax.XAxisLocation = 'top';

%%% Y ticks: one per subject, at its first recording, so the label heads the block the way
%%% the separator rule just above it does.
%%% The label font shrinks to whatever the block pitch allows, and only once even 6 pt will
%%% not fit are labels thinned out: overprinted labels are worse than fewer labels, because
%%% they are still read as if they pointed somewhere.
drawnow;
axPix = getpixelposition(ax, true);
%%% 0.9: the southoutside legend is added below and shrinks the axes a little.
pitch = 0.9 * max(1, axPix(4)) / numel(firstIdx);   % px of axis height per participant
fs    = min(10, floor(pitch / 1.35));         % ~1.35 px of line height per pt
if fs >= 6
    keep = 1:numel(firstIdx);
else
    fs   = 6;
    keep = 1:ceil(6 * 1.35 / pitch):numel(firstIdx);
end
ax.YAxis.FontSize = fs;
ax.YTick          = firstIdx(keep);
ax.YTickLabel     = subj(firstIdx(keep));

%%% Grid lines between cells, a light rule where a new participant starts, and a heavier
%%% one splitting the inputs off from the outputs.
hold(ax, 'on');
for x = 1.5:1:numel(colLabels)-0.5, xline(ax, x, 'Color', [1 1 1], 'LineWidth', 0.5); end
for k = 2:numel(firstIdx)
    yline(ax, firstIdx(k) - 0.5, 'Color', [0.6 0.6 0.6 0.4], 'LineWidth', 0.75);
end
if nInputCols > 0 && nInputCols < numel(colLabels)
    xline(ax, nInputCols + 0.5, 'Color', [0.15 0.15 0.15], 'LineWidth', 2);
end
hold(ax, 'off');

nMissIn  = nnz(S(:, 1:nInputCols) == 0);
nMissOut = nnz(S(:, nInputCols+1:end) == 0);
title(ax, sprintf('Prep pipeline status  -  %d recordings, %d missing input(s), %d missing output(s)', ...
    nRow, nMissIn, nMissOut), 'Interpreter', 'none');

%%% Legend via dummy patches
labels = {'missing', 'present', 'not checked'};
cidx   = [2 3 1];
h = gobjects(1, numel(labels));
for k = 1:numel(labels)
    h(k) = patch(ax, NaN, NaN, cmap(cidx(k), :), 'EdgeColor', 'k');
end
legend(ax, h, labels, 'Location', 'southoutside', 'Orientation', 'horizontal', 'Box', 'off');

if ~isempty(savefile)
    d = fileparts(savefile);
    if ~isempty(d) && ~exist(d, 'dir'), mkdir(d); end
    try
        exportgraphics(fig, savefile, 'Resolution', 150);
    catch
        saveas(fig, savefile);
    end
    fprintf('\nHeatmap written to %s\n', savefile);
end
end

% -------------------------------------------------------------------------
function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
