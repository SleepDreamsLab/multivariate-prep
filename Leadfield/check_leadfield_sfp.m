function report = check_leadfield_sfp(BIDS, leadfieldpath, opts)
% CHECK_LEADFIELD_SFP  Verify each head model's electrode geometry against the right .sfp.
%
%   report = check_leadfield_sfp(BIDS, leadfieldpath)
%   report = check_leadfield_sfp(BIDS, leadfieldpath, 'subjectfilter', {'sub-drop0061'})
%
%   The only thing that can actually be wrong about a head model is the electrode
%   geometry that went into it, so that is the only thing this checks: for every
%   <sub>/<ses> folder under leadfieldpath it loads the electrode positions Brainstorm
%   actually used (<sesFolder>/channel.mat's Channel array - what channel_align_auto
%   fitted and process_headmodel gained), parses every .sfp this subject has, and asks
%   which one that geometry matches. Filenames, BIDS entities and Brainstorm's own import
%   log are never used to decide a match; they cannot tell you what was actually gained
%   against, only what someone intended to import.
%
%   The match itself is a scale-free rigid (Procrustes, no reflection) fit between the two
%   labelled point sets, matched electrode-by-electrode by name: the same net placement
%   reduces to d ~ 0 once rotation, translation and any import unit rescaling are factored
%   out, while a different session's digitisation - same subject, same head, different
%   night - typically sits one to two orders of magnitude higher (validated by hand on
%   sub-drop0007/ses-t1: d = 0.0000 against ses-t2's .sfp - the file its own channel.mat
%   History says was imported - versus d = 0.0017 against ses-t3's, on 257 matched
%   channels). A session is flagged when some OTHER session's .sfp fits its geometry
%   better than its own does.
%
%   It also flags two problems that land on the same file bidsfun_gedai actually reads
%   (<sub>/<ses>/headmodel_surf_openmeeg.mat, always this exact name - see
%   gedai.loadrefcov), independent of the geometry check above:
%     - a stray headmodel_surf_openmeeg_NN.mat next to it. Brainstorm never overwrites, so
%       whenever a recompute's "delete the stale head model first" step did not run
%       (bidsfun_build_leadfield), two head models are left in the folder; bidsfun_gedai
%       always reads the unnumbered one, which is then not necessarily the one most
%       recently computed.
%     - headmodel_surf_openmeeg.mat older than channel.mat, i.e. it was computed before
%       the channel file's most recent import and so belongs to whichever geometry
%       preceded the one checked above.
%
%   Nothing under leadfieldpath is written; the only output is the CSV named by
%   opts.csvout.
%
% INPUTS:
%   BIDS           bids-matlab layout struct (e.g. BIDS_DROP{1} in PrepPipelineDROP.m) -
%                  needed only so gedai.matchSfpFile can locate every .sfp candidate.
%   leadfieldpath  Brainstorm protocol's data/ folder - the leadfieldpath variable in
%                  PrepPipelineDROP.m
%
% OPTIONAL NAME-VALUE:
%   subjectfilter      cell array of subject folder names (e.g. {'sub-drop0061'}); {} = all
%   sessionfilter       cell array of session folder names (e.g. {'ses-t1'});     {} = all
%   minmatchedchannels  electrodes a candidate .sfp must share (by name) with channel.mat
%                       before its fit is trusted; below this the comparison is skipped as
%                       inconclusive rather than risking a fit on too few points (default 20)
%   csvout              path to write the full per-session table as CSV
%                       (default <leadfieldpath>/leadfield_sfp_check.csv)
%
% OUTPUT:
%   report  table, one row per <sub>/<ses> folder found under leadfieldpath, including
%           hasSfp (does gedai.matchSfpFile resolve an actual .sfp file for THIS session -
%           independent of hasChannel/hasHeadModel; false means the montage genuinely does
%           not exist yet, not that anything was built wrong) and hasChannel/hasHeadModel
%           (do channel.mat / headmodel_surf_openmeeg.mat exist in this folder). A summary
%           of the flagged rows is also printed.

arguments
    BIDS
    leadfieldpath          char
    opts.subjectfilter     cell = {}
    opts.sessionfilter     cell = {}
    opts.minmatchedchannels (1,1) double {mustBePositive} = 20
    opts.csvout             char = ''
end

if ~isfolder(leadfieldpath)
    error('check_leadfield_sfp:noLeadfieldPath', 'leadfieldpath does not exist: %s', leadfieldpath);
end
if isempty(opts.csvout)
    opts.csvout = fullfile(leadfieldpath, 'leadfield_sfp_check.csv');
end

fprintf('\n=== Running check_leadfield_sfp ===\n');

subDirs = dir(fullfile(leadfieldpath, 'sub-*'));
subDirs = subDirs([subDirs.isdir]);

rows = {};
for iSub = 1:numel(subDirs)
    subName = subDirs(iSub).name;
    if ~isempty(opts.subjectfilter) && ~any(strcmpi(subName, opts.subjectfilter))
        continue
    end

    sesDirs = dir(fullfile(subDirs(iSub).folder, subName, 'ses-*'));
    sesDirs = sesDirs([sesDirs.isdir]);
    candSessions = {sesDirs.name};

    %%% Resolve and parse every candidate .sfp for this subject once, up front, so the
    %%% N sessions below share the same N parsed candidates instead of re-reading and
    %%% re-parsing the same files repeatedly.
    cand = buildCandidates(BIDS, subName, candSessions);

    for iSes = 1:numel(sesDirs)
        sesName = sesDirs(iSes).name;
        if ~isempty(opts.sessionfilter) && ~any(strcmpi(sesName, opts.sessionfilter))
            continue
        end
        fprintf('%s / %s ...\n', subName, sesName);
        sesFolder = fullfile(sesDirs(iSes).folder, sesName);
        rows(end+1, :) = checkOneSession(subName, sesName, sesFolder, cand, opts.minmatchedchannels); %#ok<AGROW>
    end
end

varNames = {'subject', 'session', 'hasSfp', 'hasChannel', 'hasHeadModel', 'nExtraHeadModels', ...
    'nMatchedChannels', 'ownSessionDistance', 'bestMatchSession', 'bestMatchDistance', ...
    'geometryMismatch', 'importedSfp', 'headmodelOlderThanChannel', 'flagged', 'notes'};
report = cell2table(rows, 'VariableNames', varNames);

nFlagged = sum(report.flagged);
fprintf('\n%d of %d session(s) flagged\n', nFlagged, height(report));
if nFlagged > 0
    disp(report(report.flagged, {'subject', 'session', 'notes'}))
end

writetable(report, opts.csvout);
fprintf('Full report -> %s\n', opts.csvout);
end

% -------------------------------------------------------------------------
function cand = buildCandidates(BIDS, subName, candSessions)
% One parsed {session, label, xyz} entry per session of subName that has a resolvable
% .sfp - via gedai.matchSfpFile, the same resolver the rest of the pipeline trusts, so
% this list can never itself be the source of a wrong pick. A session with no resolvable
% .sfp (the montage genuinely does not exist) is simply absent from the list; that is
% itself informative when nothing later matches a session that IS absent here.
cand = struct('session', {}, 'sfpPath', {}, 'label', {}, 'xyz', {});
for c = 1:numel(candSessions)
    ses = candSessions{c};
    try
        sfpPath = gedai.matchSfpFile(BIDS.pth, erase(subName, 'sub-'), ses);
    catch
        sfpPath = '';
    end
    if isempty(sfpPath) || ~isfile(sfpPath)
        continue
    end
    [label, xyz] = readSfp(sfpPath);
    if isempty(label)
        continue
    end
    cand(end+1) = struct('session', ses, 'sfpPath', sfpPath, 'label', {label}, 'xyz', {xyz}); %#ok<AGROW>
end
end

% -------------------------------------------------------------------------
function [label, xyz] = readSfp(sfpPath)
% .sfp is whitespace-delimited "Label X Y Z" per line, no header - fiducials (FidNz,
% FidT9, FidT10), electrodes (E1..) and the reference (Cz) all in the same format.
label = {}; xyz = zeros(0, 3);
fid = fopen(sfpPath, 'r');
if fid < 0
    return
end
C = textscan(fid, '%s %f %f %f');
fclose(fid);
if numel(C{1}) < 3
    return
end
label = C{1};
xyz   = [C{2}, C{3}, C{4}];
end

% -------------------------------------------------------------------------
function [d, nMatched] = shapeDistance(chanName, chanLoc, label, xyz, minMatched)
% Scale-free rigid (no-reflection) Procrustes distance between the channel positions
% Brainstorm actually used (chanName/chanLoc) and one candidate .sfp (label/xyz), matched
% electrode-by-electrode by name. d is bounded in [0,1]; ~0 means the same net placement,
% once rotation, translation and any uniform import rescaling are optimally factored out -
% see the module help for the validated same-session-vs-different-session gap.
[tf, iLoc] = ismember(lower(label), lower(chanName));
matched = find(tf);
nMatched = numel(matched);
if nMatched < minMatched
    d = NaN;
    return
end
X = chanLoc(iLoc(matched), :);
Y = xyz(matched, :);
d = procrustes(X, Y, 'Reflection', false);
end

% -------------------------------------------------------------------------
function row = checkOneSession(subName, sesName, sesFolder, cand, minMatched)
% One row of the report table for a single <sub>/<ses> folder. Kept as its own function
% so every exit path still returns a row of the right width - a session that turns out to
% have no channel.mat must still line up with the rest of the table.

chanFile = fullfile(sesFolder, 'channel.mat');
hmFile   = fullfile(sesFolder, 'headmodel_surf_openmeeg.mat');
hasSfp     = ~isempty(cand) && any(strcmpi({cand.session}, sesName));
hasChannel = isfile(chanFile);
hasHM      = isfile(hmFile);
nExtraHM   = numel(dir(fullfile(sesFolder, 'headmodel_surf_openmeeg_*.mat')));

importedSfp = '';
hmOlderThanChan = false;
nMatchedChannels = 0;
ownDist = NaN; bestSes = ''; bestDist = NaN;
geomMismatch = false;
notes = {};

if ~hasChannel
    notes{end+1} = 'no channel.mat';
else
    S = load(chanFile, 'Channel', 'History');

    %%% Informational only - what Brainstorm's own log says it imported. Never used to
    %%% decide the flag below; that decision is geometry-only (see module help).
    if isfield(S, 'History') && ~isempty(S.History)
        importRows = find(strcmpi(S.History(:, 2), 'import'));
        if ~isempty(importRows)
            msg = S.History{importRows(end), 3};
            tok = regexp(msg, 'Import from:\s*(.*?)\s*\(Format', 'tokens', 'once');
            if ~isempty(tok)
                importedSfp = strrep(tok{1}, '\', '/');
            end
        end
    end

    %%% The geometry check: which candidate .sfp does this channel.mat actually match?
    if isfield(S, 'Channel') && ~isempty(S.Channel)
        chanName = {S.Channel.Name};
        chanLoc  = cell2mat({S.Channel.Loc})';   % 3x1 per channel -> 3xN -> Nx3

        if isempty(cand)
            notes{end+1} = 'no candidate .sfp for this subject to compare against';
        else
            dists = nan(1, numel(cand));
            for c = 1:numel(cand)
                [dists(c), nm] = shapeDistance(chanName, chanLoc, cand(c).label, cand(c).xyz, minMatched);
                if strcmpi(cand(c).session, sesName)
                    nMatchedChannels = nm;
                end
            end
            if all(isnan(dists))
                notes{end+1} = 'no candidate .sfp shared enough channels to compare';
            else
                [bestDist, iBest] = min(dists);
                bestSes = cand(iBest).session;
                iOwn = find(strcmpi({cand.session}, sesName), 1);
                if ~isempty(iOwn)
                    ownDist = dists(iOwn);
                end
                geomMismatch = ~strcmpi(bestSes, sesName);
                if geomMismatch
                    if isnan(ownDist)
                        notes{end+1} = sprintf(['channel geometry matches %s''s .sfp (d=%.4g); ' ...
                            'no .sfp exists for %s itself'], bestSes, bestDist, sesName);
                    else
                        notes{end+1} = sprintf(['channel geometry matches %s''s .sfp (d=%.4g) better ' ...
                            'than its own session''s (d=%.4g)'], bestSes, bestDist, ownDist);
                    end
                end
            end
        end
    end

    if hasHM
        dc = dir(chanFile); dh = dir(hmFile);
        if dh.datenum < dc.datenum
            hmOlderThanChan = true;
            notes{end+1} = 'headmodel predates the current channel import';
        end
    end
end

if nExtraHM > 0
    notes{end+1} = sprintf('%d stale headmodel_surf_openmeeg_NN.mat alongside it', nExtraHM);
end

flagged = geomMismatch || hmOlderThanChan || (nExtraHM > 0);

row = {subName, sesName, hasSfp, hasChannel, hasHM, nExtraHM, nMatchedChannels, ownDist, bestSes, bestDist, ...
    geomMismatch, importedSfp, hmOlderThanChan, flagged, strjoin(notes, '; ')};
end
