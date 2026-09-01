function build_leadfield_bids(bids, opts)
% BUILD_LEADFIELD_BIDS  Scale-only fsaverage leadfields from BIDS layout + .sfp files.
%
% USAGE:
%   build_leadfield_bids(bids)
%   build_leadfield_bids(bids, SubjectFilter={'sub-hpmam003','sub-hpmam004'})
%   build_leadfield_bids(bids, DoQC=false, nScalp=642)
%
% INPUTS:
%   bids  — bids-matlab layout struct (e.g. BIDS_PM{1}). The study type
%           ('PM', 'DROP', or unknown) is inferred from bids.pth, and the
%           sfp search root is derived from it accordingly.
%
% OPTIONAL NAME-VALUE (opts):
%   SubjectFilter    cell array of chars listing subject IDs to process; {} = all.
%   SfpResolver      function handle (bidsPath, subjectName, sessName, studyId)
%                    -> sfpPath. Defaults to the built-in bidsToSfp, which
%                    derives the search root from bids.pth per study:
%                      'PM'   — fileparts(bids.pth)/Data_collection/... [VERIFY]
%                      'DROP' — bids.pth/sourcedata/gps/.../solved/*.sfp,
%                               then bids.pth/rawdata/<sub>/<ses>/eeg/*.sfp;
%                               candidates naming a different ses-* are
%                               rejected, so a session with no file is
%                               skipped rather than given another's coords
%                      other  — returns bids.pth/<sub>/<ses>/eeg/ so skip
%                               messages show where to look
%                    Supply your own to support a custom folder layout, e.g.:
%                      opts.SfpResolver = @(p,sub,ses,~) fullfile(p, sub, ses, 'coords.sfp')
%   ProtocolName     Brainstorm protocol name         (default 'GEDAI_Leadfield')
%   TemplateName     fsaverage template name          (default 'FsAverage_2020')
%   SfpFormat        channel-file format string       (default 'EGI')
%   nScalp/nOuter/nInner  BEM mesh vertex counts     (default 1922 each)
%   WarpTolerance    fraction of head points dropped as outliers before
%                    scale warp; 0 = keep all; must be in [0,1)  (default 0)
%   iWarpRefSession  1-based index of session driving the scale warp (default 1)
%   ForceReprocess   reprocess subjects that already have BEM surfaces; when
%                    false (default) such subjects are skipped; when true the
%                    anatomy is reset to the template before re-warping so the
%                    scale is not applied twice                   (default false)
%   DoQC             save registration PNG figures               (default true)
%   QCDir            output folder for QC images      (default <pwd>/QC_registration)

arguments
    bids             struct
    opts.SubjectFilter   (1,:) cell    = {}
    opts.SfpResolver     (1,1) function_handle = @bidsToSfp
    opts.ProtocolName    (1,1) string  = "GEDAI_Leadfield"
    opts.TemplateName    (1,1) string  = "FsAverage_2020"
    opts.SfpFormat       (1,1) string  = "EGI"
    opts.nScalp          (1,1) double {mustBePositive,mustBeInteger} = 1922
    opts.nOuter          (1,1) double {mustBePositive,mustBeInteger} = 1922
    opts.nInner          (1,1) double {mustBePositive,mustBeInteger} = 1922
    opts.WarpTolerance   (1,1) double {mustBeNonnegative}            = 0
    opts.iWarpRefSession (1,1) double {mustBePositive,mustBeInteger} = 1
    opts.ForceReprocess  (1,1) logical = false
    opts.DoQC            (1,1) logical = true
    opts.QCDir           (1,1) string  = ""
end
if opts.WarpTolerance >= 1
    error('WarpTolerance must be in [0,1); got %.4g', opts.WarpTolerance);
end
if opts.QCDir == ""
    opts.QCDir = string(fullfile(pwd, 'QC_registration'));
end

% Auto-detect study from the BIDS root path; used by the default bidsToSfp.
if contains(bids.pth, 'data-drop', 'IgnoreCase', true)
    studyId = 'DROP';
elseif ~isempty(regexp(bids.pth, '[/\\]PM([/\\]|$)', 'once'))
    studyId = 'PM';
else
    studyId = '';
end

openmeegOpt = struct( ...
    'BemFiles', {{}}, 'BemNames', {{'Scalp','Skull','Brain'}}, ...
    'BemCond', [1, 0.0125, 1], 'BemSelect', [1, 1, 1], ...
    'isAdjoint', 0, 'isAdaptative', 1, 'isSplit', 0, 'SplitLength', 4000);

%% Protocol setup
if ~brainstorm('status'); brainstorm nogui; end
iProtocol = bst_get('Protocol', char(opts.ProtocolName));
if isempty(iProtocol)
    gui_brainstorm('CreateProtocol', char(opts.ProtocolName), 1, 0);
else
    gui_brainstorm('SetCurrentProtocol', iProtocol);
end
bst_report('Start');
if opts.DoQC && ~exist(char(opts.QCDir), 'dir')
    mkdir(char(opts.QCDir));
end

%% Set FsAverage as default anatomy for the protocol.
sTemplates = bst_get('AnatomyDefaults');
iTemplate  = find(strcmpi(char(opts.TemplateName), {sTemplates.Name}), 1);
if isempty(iTemplate)
    error('Template "%s" not found. Available: %s', opts.TemplateName, strjoin({sTemplates.Name}, ', '));
end
db_set_template(0, sTemplates(iTemplate), 0);
db_save();

%% Group BIDS entries by participant; apply subject filter if set.
names  = {bids.subjects.name};
uNames = unique(names, 'stable');
if ~isempty(opts.SubjectFilter)
    uNames = uNames(contains(uNames, opts.SubjectFilter));
end

for p = 1:numel(uNames)
    subjectName = uNames{p};
    iEntries    = find(strcmp(names, subjectName));   % indices of all sessions for this subject

    % Resolve a .sfp for each session; skip sessions with none.
    sessName = {}; sfpPath = {};
    for e = 1:numel(iEntries)
        sess = bids.subjects(iEntries(e)).session;
        if ~isempty(sess)
            try
                sfp = opts.SfpResolver(bids.pth, subjectName, sess, studyId);
            catch ME
                fprintf('[skip] %s / %s: SfpResolver error — %s\n', subjectName, sess, ME.message);
                sfp = '';
            end
            if ~isempty(sfp) && isfile(sfp)
                sessName{end+1} = sess;   %#ok<AGROW>
                sfpPath{end+1}  = sfp;    %#ok<AGROW>
            elseif isempty(sfp)
                fprintf('[skip] %s / %s: no .sfp found\n', subjectName, sess);
            else
                fprintf('[skip] %s / %s: no .sfp (looked for: %s)\n', subjectName, sess, sfp);
            end
        end
    end
    if isempty(sfpPath)
        fprintf('[skip subject] %s: no .sfp for any session\n', subjectName);
        continue;
    end
    nSess = numel(sfpPath);

    % Detect whether a prior run already completed warp + BEM for this subject.
    % BEM surface filenames contain 'bem', making them a reliable proxy.
    [sSubjectPre, iSubjectPre] = bst_get('Subject', subjectName, 0);
    hasBEM = ~isempty(iSubjectPre) && iSubjectPre > 0 && ...
        ~isempty(sSubjectPre) && ~isempty(sSubjectPre.Surface) && ...
        any(contains({sSubjectPre.Surface.FileName}, 'bem'));

    if hasBEM && ~opts.ForceReprocess
        fprintf('[skip] %s: BEM surfaces already exist (set ForceReprocess=true to redo)\n', subjectName);
        continue;
    end

    % Create subject if absent; guards re-runs without triggering dialogs.
    if isempty(iSubjectPre) || iSubjectPre == 0
        db_add_subject(subjectName, [], 1, 0);
    elseif hasBEM
        % Force-reprocessing an already-warped subject: reset the anatomy back
        % to the template so bst_warp_prepare starts from unscaled surfaces.
        fprintf('[reset] %s: resetting anatomy to template before re-warp\n', subjectName);
        db_set_template(iSubjectPre, sTemplates(iTemplate), 0);   % 0 = non-interactive: skip the java_dialog confirm that hangs in nogui
        db_save();
    end

    %% Per-session: import .sfp and refine electrode registration (rigid ICP).
    sessStudy = zeros(nSess, 1);
    sessChan  = cell(nSess, 1);
    for s = 1:nSess
        iStudy = db_add_condition(subjectName, sessName{s});
        % import_channel(iStudies, File, Format, ChannelReplace, ChannelAlign, isSave, isFixUnits, isApplyVox2ras)
        import_channel(iStudy, sfpPath{s}, char(opts.SfpFormat), 2, 0, 1, 1, 0);

        sStudy      = bst_get('Study', iStudy);
        channelFile = sStudy.Channel.FileName;

        % Rigid ICP refinement. Electrodes are NOT added as head points:
        % bst_warp merges them internally; duplicating them breaks the
        % warp's outlier-removal indexing.
        channelMat = channel_align_auto(channelFile, [], 0, 0);
        if ~isempty(channelMat)
            bst_save(file_fullpath(channelFile), channelMat, 'v7');
        end

        sessStudy(s) = iStudy;
        sessChan{s}  = channelFile;

        % QC 1: pre-warp electrode fit (green = close to scalp, red = far).
        if opts.DoQC
            hFig = channel_align_manual(channelFile, 'EEG', 0);
            for view = {'left', 'front', 'top'}
                figure_3d('SetStandardView', hFig, view{1});
                drawnow;
                qcFile = char(fullfile(opts.QCDir, sprintf('pre_warp_%s_%s_%s.png', subjectName, sessName{s}, view{1})));
                saveas(hFig, qcFile);
                fprintf('[QC pre-warp] %s\n', qcFile);
            end
            close(hFig);
        end
    end

    %% Scale-only warp of fsaverage to the reference session.
    iRef    = min(opts.iWarpRefSession, nSess);
    warpOpt = struct('isScaleOnly', 1, 'tolerance', opts.WarpTolerance, ...
                     'isSurfaceOnly', 0, 'isInterp', 1, 'isInteractive', 0);
    bst_warp_prepare(sessChan{iRef}, warpOpt);
    db_save();

    % QC 2: post-warp electrode fit per session on the shared scaled scalp.
    if opts.DoQC
        sSubjQC   = bst_get('Subject', subjectName);
        scalpFile = sSubjQC.Surface(sSubjQC.iScalp).FileName;
        for s = 1:nSess
            hFig = view_headpoints(sessChan{s}, scalpFile);
            for view = {'left','front','top'}
                figure_3d('SetStandardView', hFig, view{1});
                drawnow;
                qcFile = char(fullfile(opts.QCDir, sprintf('post_warp_%s_%s_%s.png', subjectName, sessName{s}, view{1})));
                saveas(hFig, qcFile);
                fprintf('[QC post-warp] %s\n', qcFile);
            end
            close(hFig);
        end
    end

    %% Generate BEM surfaces on the scaled anatomy.
    bst_process('CallProcess', 'process_generate_bem', [], [], ...
        'subjectname', subjectName, ...
        'nscalp', opts.nScalp, 'nouter', opts.nOuter, 'ninner', opts.nInner, ...
        'thickness', 4, 'method', 'brainstorm');
    db_save();

    %% Per-session: compute OpenMEEG head model and report Gain matrix size.
    for s = 1:nSess
        % The BEM surfaces were just rebuilt, so any head model still sitting in
        % this folder is stale. Brainstorm never overwrites — bst_headmodeler
        % passes the output name through file_unique — so leaving the old file
        % in place yields headmodel_surf_openmeeg_02.mat, _03.mat, ... Delete
        % first so the recomputed model lands on the canonical filename.
        [sStudyOld, iStudyOld] = bst_get('ChannelFile', sessChan{s});
        if ~isempty(iStudyOld) && ~isempty(sStudyOld) && ~isempty(sStudyOld.HeadModel)
            for h = 1:numel(sStudyOld.HeadModel)
                oldFile = file_fullpath(sStudyOld.HeadModel(h).FileName);
                if exist(oldFile, 'file')
                    file_delete(oldFile, 1);
                end
            end
            fprintf('[clean] %s / %s: removed %d stale head model(s)\n', ...
                subjectName, sessName{s}, numel(sStudyOld.HeadModel));
            sStudyOld.HeadModel  = repmat(db_template('HeadModel'), 0, 1);
            sStudyOld.iHeadModel = [];
            bst_set('Study', iStudyOld, sStudyOld);
            db_save();
        end

        nRep = reportRowCount();
        bst_process('CallProcess', 'process_headmodel', [], [], ...
            'channelfile', sessChan{s}, ...
            'sourcespace', 1, 'meg', 1, 'eeg', 3, ...
            'openmeeg', openmeegOpt);

        % Re-fetch study from the channel file rather than the stored index,
        % which can go stale after BEM/warp steps reorganise the database.
        [sStudyHM, ~] = bst_get('ChannelFile', sessChan{s});
        if isempty(sStudyHM) || isempty(sStudyHM.iHeadModel)
            fprintf('[warn] %s / %s: head model missing after computation\n', subjectName, sessName{s});
            % process_headmodel reports its reason to the Brainstorm report and
            % returns silently; echo it so the console says why it gave up.
            msgs = reportMessagesSince(nRep);
            for m = 1:numel(msgs)
                fprintf('        %s\n', msgs{m});
            end
            if isempty(msgs)
                fprintf('        (no error recorded in the Brainstorm report)\n');
            end
            continue;
        end
        headModelFile = sStudyHM.HeadModel(sStudyHM.iHeadModel).FileName;
        hm            = in_bst_headmodel(headModelFile);
        fprintf('[%s / %s] Gain %d x %d (x3). %s\n', ...
            subjectName, sessName{s}, size(hm.Gain,1), size(hm.GridLoc,1), headModelFile);
    end
end

reportFile = bst_report('Save', []);
bst_report('Export', reportFile, fullfile(pwd, 'leadfield_bids_report.html'));
fprintf('Done. Report: %s\n', fullfile(pwd, 'leadfield_bids_report.html'));
end


% -------------------------------------------------------------------------
function n = reportRowCount()
% Current number of entries in the running Brainstorm process report.
global GlobalData
n = 0;
try
    n = size(GlobalData.ProcessReports.Reports, 1);
catch
end
end

% -------------------------------------------------------------------------
function msgs = reportMessagesSince(nBefore)
% Error/warning messages appended to the Brainstorm process report after row
% nBefore. Column 1 is the entry type, column 4 the message text.
global GlobalData
msgs = {};
try
    R = GlobalData.ProcessReports.Reports;
catch
    return;
end
if isempty(R) || size(R, 2) < 4
    return;
end
for i = (nBefore + 1):size(R, 1)
    if any(strcmpi(R{i,1}, {'error', 'warning'}))
        txt = strtrim(strrep(char(R{i,4}), char(10), ' | '));
        msgs{end+1} = sprintf('%s: %s', upper(char(R{i,1})), txt);   %#ok<AGROW>
    end
end
end

% -------------------------------------------------------------------------
function sfp = bidsToSfp(bidsPath, subjectName, sessName, studyId)
% Dispatch .sfp lookup based on studyId. [VERIFY] Edit the sub-functions
% below if your naming or folder layout differs.
%
%   'PM'   — fileparts(bidsPath)/Data_collection/<SubColl>/SA_stim/<ses>/GPS/*.sfp
%   'DROP' — bidsPath/sourcedata/gps/<sub>/solved/*domesolved*.sfp,
%             then bidsPath/rawdata/<sub>/<ses>/eeg/*.sfp
%   other  — returns bidsPath/<sub>/<ses>/eeg/ so skip messages show where to look

switch upper(studyId)
    case 'PM'
        sfp = sfpFromPM(bidsPath, subjectName, sessName);
    case 'DROP'
        sfp = sfpFromDrop(bidsPath, subjectName, sessName);
    otherwise
        sfp = fullfile(bidsPath, subjectName, sessName, 'eeg');
end
end

% -------------------------------------------------------------------------
function sfp = sfpFromPM(bidsPath, subjectName, sessName)
% PM: Data_collection is assumed to be a sibling of the BIDS folder. [VERIFY]
    collRoot = fullfile(fileparts(bidsPath), 'Data_collection');
    lab = erase(subjectName, 'sub-');
    num = regexp(lab, '\d+$', 'match', 'once');
    if isempty(num), sfp = ''; return; end
    subColl = sprintf('H%s_PM_AM', num);   % [VERIFY skeleton]
    ses     = erase(sessName, 'ses-');
    gpsDir  = fullfile(collRoot, subColl, 'SA_stim', ses, 'GPS');
    sfp     = fullfile(gpsDir, sprintf('%s_GPS_%s_coordinates.sfp', subColl, ses));
    if ~isfile(sfp)
        d = dir(fullfile(gpsDir, '*coordinates*.sfp'));
        if isempty(d), d = dir(fullfile(gpsDir, '*.sfp')); end
        if ~isempty(d), sfp = fullfile(d(1).folder, d(1).name); end
    end
end

% -------------------------------------------------------------------------
function sfp = sfpFromDrop(bidsPath, subjectName, sessName)
% DROP: primary path is sourcedata/gps/<sub>/solved/*domesolved*.sfp;
% falls back to rawdata/<sub>/<ses>/eeg/*.sfp.
    solvedDir = fullfile(bidsPath, '..', 'sourcedata', 'gps', subjectName, 'solved');
    sessPat = ['*' sessName '*domesolved*.sfp'];
    d = dir(fullfile(solvedDir, sessPat));
    if isempty(d), d = dropForeignSessions(dir(fullfile(solvedDir, '*.sfp')), sessName); end
    if ~isempty(d)
        sfp = fullfile(d(1).folder, d(1).name);
        return;
    end
    eegDir = fullfile(bidsPath, 'rawdata', subjectName, sessName, 'eeg');
    d = dropForeignSessions(dir(fullfile(eegDir, '*.sfp')), sessName);
    if ~isempty(d)
        sfp = fullfile(d(1).folder, d(1).name);
    else
        sfp = '';
    end
end

% -------------------------------------------------------------------------
function d = dropForeignSessions(d, sessName)
% Remove candidates whose filename carries a ses-* token for a *different*
% session. Files with no session token are kept: some layouts store a single
% .sfp per subject. Without this filter the wildcard fallback silently hands
% a session whichever .sfp dir() happened to list first, which can be another
% session's electrode positions.
    keep = true(size(d));
    for i = 1:numel(d)
        tok = regexp(d(i).name, 'ses-[A-Za-z0-9]+', 'match');
        if ~isempty(tok) && ~any(strcmpi(tok, sessName))
            keep(i) = false;
        end
    end
    d = d(keep);
end
