function build_leadfield_bids(BIDS, CollectionRoot)
% BUILD_LEADFIELD_BIDS  Data-free scaled-fsaverage leadfields, driven by a
% bids-matlab layout, with electrode positions pulled from a separate
% Data_collection tree of .sfp files.
%
% USAGE:
%   BIDS = BIDS_PM{1};
%   build_leadfield_bids(BIDS, 'W:\PM\Data_collection');
%
% The BIDS struct only enumerates subjects/sessions and provides the
% Brainstorm naming. Geometry comes entirely from the mapped .sfp; the
% _eeg.dat recordings are never read.
%
% Per participant (grouped from BIDS.subjects by .name):
%   per session: condition + import .sfp (fiducial SCS) + refine (rigid ICP)
%   once:        SCALE-ONLY warp of fsaverage to a reference session + BEM
%   per session: OpenMEEG BEM head model via 'channelfile' -> Gain
%
% [VERIFY] tags: (1) BIDS<->collection subject mapping in bidsToSfp,
%                (2) channel_align_auto saves back,
%                (3) bst_warp_prepare runs non-interactively.

%% ===================== CONFIG =====================
ProtocolName    = 'GEDAI_Leadfield';
TemplateName    = 'FsAverage_2020';
SfpFormat       = 'EGI';          % EGI GPS .sfp
nScalp = 1922; nOuter = 1922; nInner = 1922;
% Outlier head points dropped before scaling. NOTE: bst_warp_prepare uses
% nRemove = ceil(WarpTolerance * nPoints), so this is a FRACTION (0-1), NOT a
% percentage. 0 keeps all points (and skips the removal block). Do not set >=1.
WarpTolerance = 0;
iWarpRefSession  = 1;             % which session (in listed order) drives the scaling
DoQC             = true;          % save registration figures automatically (set false to skip)
QCDir            = fullfile(pwd, 'QC_registration');  % folder for saved QC images

OpenmeegOpt = struct( ...
    'BemFiles', {{}}, 'BemNames', {{'Scalp','Skull','Brain'}}, ...
    'BemCond', [1, 0.0125, 1], 'BemSelect', [1, 1, 1], ...
    'isAdjoint', 0, 'isAdaptative', 1, 'isSplit', 0, 'SplitLength', 4000);

%% ===================== PROTOCOL =====================
if ~brainstorm('status'); brainstorm nogui; end
iProtocol = bst_get('Protocol', ProtocolName);
if isempty(iProtocol)
    gui_brainstorm('CreateProtocol', ProtocolName, 1, 0);
else
    gui_brainstorm('SetCurrentProtocol', iProtocol);
end
bst_report('Start');
if DoQC && ~exist(QCDir,'dir'); mkdir(QCDir); end

%% ===================== DEFAULT ANATOMY = FsAverage =====================
sTemplates = bst_get('AnatomyDefaults');
iTemplate  = find(strcmpi(TemplateName, {sTemplates.Name}), 1);
if isempty(iTemplate)
    error('Template "%s" not found. Available: %s', TemplateName, strjoin({sTemplates.Name}, ', '));
end
db_set_template(0, sTemplates(iTemplate), 0);
db_save();

%% ===================== GROUP BIDS subjects x sessions BY PARTICIPANT ====
names = {BIDS.subjects.name};
[uNames, ~, grp] = unique(names, 'stable');

for p = 3%1:numel(uNames)
    SubjectName = uNames{p};            % e.g. 'sub-hpmam003' -> Brainstorm subject
    iEntries    = find(grp == p);       % one entry per session

    % --- Resolve a .sfp for each session (skip sessions with none) ---
    SessName = {}; SfpPath = {};
    for e = 1:numel(iEntries)
        sess = BIDS.subjects(iEntries(e)).session;     % 'ses-S1'
        sfp  = bidsToSfp(CollectionRoot, SubjectName, sess);
        if ~isempty(sfp) && isfile(sfp)
            SessName{end+1} = sess;     %#ok<AGROW>
            SfpPath{end+1}  = sfp;      %#ok<AGROW>
        else
            fprintf('[skip] %s / %s: no .sfp found (looked for: %s)\n', SubjectName, sess, sfp);
        end
    end
    if isempty(SfpPath)
        fprintf('[skip subject] %s: no .sfp for any session\n', SubjectName);
        continue;
    end
    nSess = numel(SfpPath);

    % Only create the subject if it doesn't already exist in the protocol
    % (guards against re-running a partial batch without dialog popups).
    [~, iSubjectCheck] = bst_get('Subject', SubjectName, 0);
    if isempty(iSubjectCheck) || iSubjectCheck == 0
        db_add_subject(SubjectName, [], 1, 0);   % start on fsaverage, own channels
    end

    % --- Import + head points + refine, per session ---
    SessStudy = zeros(nSess,1);
    SessChan  = cell(nSess,1);
    for s = 1:nSess
        iStudy = db_add_condition(SubjectName, SessName{s});

        % import_channel(iStudies, File, Format, ChannelReplace, ChannelAlign, isSave, isFixUnits, isApplyVox2ras)
        import_channel(iStudy, SfpPath{s}, SfpFormat, 2, 0, 1, 1, 0);  % 2=silent overwrite

        sStudy      = bst_get('Study', iStudy);
        ChannelFile = sStudy.Channel.FileName;

        % REFINE registration (rigid ICP; layout preserved). Do NOT add the
        % electrodes as head points: for EEG, refine and warp augment the
        % registration with the EEG electrode locations automatically. Adding
        % them duplicates the set and breaks the warp's de-duplication
        % ("Remove N points" / index-exceeds error).
        % [VERIFY] channel_align_auto saves back; we re-save to be safe.
        ChannelMat = channel_align_auto(ChannelFile, [], 0, 0);
        if ~isempty(ChannelMat)
            bst_save(file_fullpath(ChannelFile), ChannelMat, 'v7');
        end

        SessStudy(s) = iStudy;
        SessChan{s}  = ChannelFile;

        % [QC 1] Pre-warp: check refine fit (isEdit=0 = view only, not editable).
        % Shows the (unscaled) scalp + electrodes colour-coded by distance.
        % Green = close to scalp, red = far. Run this for session 1 of each
        % subject at minimum; it catches fiducial errors before the warp.
        if DoQC
            hFigQC1 = channel_align_manual(ChannelFile, 'EEG', 0);
            qcFile1 = fullfile(QCDir, sprintf('pre_warp_%s_%s.png', SubjectName, SessName{s}));
            saveas(hFigQC1, qcFile1);
            close(hFigQC1);
            fprintf('[QC pre-warp] saved %s\n', qcFile1);
        end
    end

    % --- SCALE-ONLY warp of fsaverage to the reference session ---
    iRef = min(iWarpRefSession, nSess);
    % [VERIFY] bst_warp_prepare applies non-interactively. If it only previews,
    %          run this one step in the GUI (with visual QC) per subject.
    WarpOpt = struct('isScaleOnly', 1, 'tolerance', WarpTolerance, ...
                     'isSurfaceOnly', 0, 'isInterp', 1, 'isInteractive', 0);
    bst_warp_prepare(SessChan{iRef}, WarpOpt);
    db_save();

    % [QC 2] Post-warp: head points on the scaled scalp.
    % The cap should now sit flush on the warped surface. If the cap still
    % floats far off the scalp, the scaling failed (bad fiducials or units).
    if DoQC
        sSubjQC   = bst_get('Subject', SubjectName);
        ScalpFile = sSubjQC.Surface(sSubjQC.iScalp).FileName;
        hFigQC2   = view_headpoints(SessChan{iRef}, ScalpFile);
        % Save all 6 standard views (keys 1-6 in the Brainstorm GUI).
        qcViews = {'left','right','front','back','top','bottom'};
        for iV = 1:numel(qcViews)
            figure_3d('SetStandardView', hFigQC2, qcViews{iV});
            drawnow;  % flush render before capture
            qcFile2 = fullfile(QCDir, sprintf('post_warp_%s_%s.png', SubjectName, qcViews{iV}));
            saveas(hFigQC2, qcFile2);
            fprintf('[QC post-warp] saved %s\n', qcFile2);
        end
        close(hFigQC2);
    end

    % --- BEM on the scaled anatomy ---
    bst_process('CallProcess', 'process_generate_bem', [], [], ...
        'subjectname', SubjectName, ...
        'nscalp', nScalp, 'nouter', nOuter, 'ninner', nInner, ...
        'thickness', 4, 'method', 'brainstorm');
    db_save();

    % --- Head model per session (no data file: use 'channelfile') ---
    for s = 1:nSess
        bst_process('CallProcess', 'process_headmodel', [], [], ...
            'channelfile', SessChan{s}, ...
            'sourcespace', 1, 'meg', 1, 'eeg', 3, ...
            'openmeeg', OpenmeegOpt);

        sStudyHM      = bst_get('Study', SessStudy(s));
        HeadModelFile = sStudyHM.HeadModel(sStudyHM.iHeadModel).FileName;
        HM            = in_bst_headmodel(HeadModelFile);

        Gain    = HM.Gain;        % [nChan x 3*nSrc], V/(A.m), unconstrained
        GridLoc = HM.GridLoc;     % [nSrc x 3]
        fprintf('[%s / %s] Gain %d x %d (x3). %s\n', ...
            SubjectName, SessName{s}, size(Gain,1), size(GridLoc,1), HeadModelFile);
    end
end

ReportFile = bst_report('Save', []);
bst_report('Export', ReportFile, fullfile(pwd, 'leadfield_bids_report.html'));
fprintf('Done. Report: %s\n', fullfile(pwd, 'leadfield_bids_report.html'));
end


% ======================================================================
function sfp = bidsToSfp(CollectionRoot, SubjectName, SessName)
% Map a BIDS (subject, session) to its .sfp in the Data_collection tree.
%
% [VERIFY] Inferred from two (non-matching) examples:
%     sub-hpmam003  ->  H003_PM_AM           (trailing number fills the skeleton)
%     ses-S1        ->  S1
%     .../H###_PM_AM/SA_stim/S#/GPS/H###_PM_AM_GPS_S#_coordinates.sfp
% Edit this single function if your naming differs.

    lab = erase(SubjectName, 'sub-');               % 'hpmam003'
    num = regexp(lab, '\d+$', 'match', 'once');     % '003'
    if isempty(num), sfp = ''; return; end
    subColl = sprintf('H%s_PM_AM', num);            % 'H003_PM_AM'  [VERIFY skeleton]

    ses = erase(SessName, 'ses-');                  % 'S1'
    gpsDir = fullfile(CollectionRoot, subColl, 'SA_stim', ses, 'GPS');

    % Preferred exact name, then glob fallback for naming variations.
    sfp = fullfile(gpsDir, sprintf('%s_GPS_%s_coordinates.sfp', subColl, ses));
    if ~isfile(sfp)
        d = dir(fullfile(gpsDir, '*coordinates*.sfp'));
        if isempty(d), d = dir(fullfile(gpsDir, '*.sfp')); end
        if ~isempty(d)
            sfp = fullfile(d(1).folder, d(1).name);
        end
    end
end