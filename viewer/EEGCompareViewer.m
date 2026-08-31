%% EEGCompareViewer
% Overlay two EEG recordings channel-by-channel for direct visual
% comparison, with epoch/stage navigation driven by a sleep scoring file.
% Requires EEGLAB (pop_interp, pop_select, ...) already on the path.
clc;
run('..\dependancies.m')

%%% EEG pathing
fileroot    = '\\vs03.herseninstituut.knaw.nl\VS03-SandC-1\data\nin\data-drop\rawdata\derivatives\prep-zc-ged\sub-drop0001\ses-t1\';
filename    = 'sub-drop0001_ses-t1_task-sleep_run-01_desc-hpzcged_eeg.set';  

%%% ICA pathing
icaroot = fileroot;
icafile = 'sub-drop0001_ses-t1_task-sleep_run-01_desc-pamica_ica';

%%% Scoring file ('' if none)
scoringFile = '\\vs03.herseninstituut.knaw.nl\VS03-SandC-1\data\nin\data-drop\derivatives\scores\final\sub-drop0001_ses-t1_task-sleep_run-01_eeg.csv';

%%% Variables
icaSource = 'pylabel';  % 'mod' = byte-identical Fortran output (loadmodout15), cross-checked against pam
                    % 'pam' = pamica-computed (.mat) decomposition alone, no Fortran binary needed
net = 'EGI256';     % net type, passed to chans1020
targetrate = 250;   % down-sample for speed
plotDecimation = 5; % only plot every Nth sample (display only; does not affect underlying data)



%% ===================== IMPORT + PREP =====================
file = fullfile(fileroot, filename);
EEG1 = fast_eeg_import(fullfile(fileroot, filename));

% Downsample
if targetrate < EEG1.srate; 
    EEG1 = pop_resample(EEG1, targetrate); 
end
EEG = EEG1;

% Load ICA file
switch icaSource
    case 'mod'
        pam = load(fullfile(icaroot, [icafile '.mat']));        
        mod = loadmodout15(fullfile(icaroot, icafile));   % W, S — byte-identical to Fortran
        EEG = icaFromMod(EEG, pam, mod, file, icafile);

        % ICLabel
        EEG = iclabel(EEG);
        EEG = selectcomps(EEG, 'ArtefactThreshold', .0, 'ManualQC', false, 'ICLabelClasses', [2:7]);
        badComps = find(EEG.reject.gcompreject);
        % EEG2 = pop_subcomp(EEG, badComps, 0);  

    case 'pam'
        pam = load(fullfile(icaroot, [icafile '.mat']));        
        EEG = icaFromPam(EEG, pam, file, icafile);

        % ICLabel
        EEG = iclabel(EEG);
        EEG = selectcomps(EEG, 'ArtefactThreshold', .0, 'ManualQC', false, 'ICLabelClasses', [2:7]);
        badComps = find(EEG.reject.gcompreject);
        % EEG2 = pop_subcomp(EEG, badComps, 0);     
        
    case 'pylabel'
        EEG = loadica(EEG, fullfile(icaroot, [icafile '.mat']))
end



% Scoring
scoringDigits = [];
if ~isempty(scoringFile)
    scoringDigits = scoreloader(scoringFile);
end

% Naming
[~, fname1] = fileparts(filename);
file = fname1;
file2 = [fname1 ' + ICA'];

% Add EOG, interpolate missing channels
chanOpts       = {false, 'add_eog', true, 'net', net, 'chanprefix', 'E'};
EEGpreChans    = EEG;                       % keep the pre-chans1020 channel set/labels
[EEG, chanmap] = chans1020(EEG, chanOpts{:});

% EOG1/EOG2 are bipolar derivations that chans1020 appends AFTER the ICA was
% built, so they sit outside icachansind and IC subtraction would leave them
% untouched -- eye components would stay fully visible in exactly the two
% traces you use to judge them. Extend the mixing matrix to cover them.
EEG = viewer.extendICAToDerived(EEG, EEGpreChans, chanOpts, [chanmap.EOG1, chanmap.EOG2]);

% View dataa
viewer.eegCompareViewer(EEG, [], scoringDigits, file, file2, plotDecimation, fieldnames(chanmap));




%% ===================== LOCAL FUNCTIONS =====================
function EEG = icaFromMod(EEG, pam, mod, file1, icafile)
% ICAFROMMOD Build EEG ICA fields from the byte-identical Fortran AMICA
% output (mod, via loadmodout15). pam supplies only P (icasphere) and
% chanlabels; pam.icaweights is used solely as an independent cross-check.
    labels = string(pam.chanlabels(:));
    P = double(pam.icasphere);                        % k x nchan, the rank projection
    W = double(mod.W(:,:,1));                         % explicit model index
    S = double(mod.S(1:mod.num_pcs, 1:mod.data_dim)); % doc: S is num_pcs x data_dim

    assert(size(P,2) == numel(labels), 'P columns vs chanlabels mismatch');
    assert(size(S,2) == size(P,1),     'amicaout rank vs P mismatch');

    % amicaout lives in the P-space, so fold P into the sphere
    mixFull   = pinv(W * S * P);                       % nchan_amica x ncomp, the mixing matrix
    [tf, idx] = ismember(labels, string({EEG.chanlocs.labels}'));

    checkICAProvenance(file1, icafile);

    if all(tf)
        % Exact path: every channel the decomposition needs is still present.
        EEG.icaweights  = W;
        EEG.icasphere   = S * P;
        EEG.icawinv     = mixFull;
        EEG.icachansind = idx(:)';

        % cross-check the two export paths (scale-invariant, so normalization is moot)
        C = abs(corr((EEG.icaweights*EEG.icasphere)', (double(pam.icaweights)*P)'));
        r = max(C, [], 2);
        fprintf('mod vs mat: median %.4f, min %.4f, n<0.99 = %d/%d\n', ...
                median(r), min(r), sum(r < 0.99), numel(r));
    else
        EEG = applyICASubset(EEG, mixFull, tf, idx, labels);
    end

    EEG = eeg_checkset(EEG, 'ica');
end

function EEG = icaFromPam(EEG, pam, file1, icafile)
% ICAFROMPAM Build EEG ICA fields from the pamica-computed (.mat)
% decomposition alone -- no Fortran binary needed. pam.icaweights already
% has the internal AMICA sphere folded in (see run-pamica.py), so icasphere
% here is just P, and icawinv is already the precomputed mixing matrix. This
% is the pamica-side (single precision) decomposition, not the byte-identical
% Fortran one that icaFromMod uses.
    labels  = string(pam.chanlabels(:));
    mixFull = double(pam.icawinv);   % nchan x ncomp mixing matrix, already pinv(icaweights*icasphere)

    assert(size(mixFull,1) == numel(labels), 'icawinv rows vs chanlabels mismatch');

    [tf, idx] = ismember(labels, string({EEG.chanlocs.labels}'));

    checkICAProvenance(file1, icafile);

    if all(tf)
        % Exact path: every channel the decomposition needs is still present.
        EEG.icaweights  = double(pam.icaweights);
        EEG.icasphere   = double(pam.icasphere);
        EEG.icawinv     = mixFull;
        EEG.icachansind = idx(:)';
    else
        EEG = applyICASubset(EEG, mixFull, tf, idx, labels);
    end

    EEG = eeg_checkset(EEG, 'ica');
end

function EEG = applyICASubset(EEG, mixFull, tf, idx, labels)
% APPLYICASUBSET Some channels the decomposition was trained on are gone --
% a different derivative dropped them as bad. The exact unmixing cannot be
% formed, because it needs every column. Estimate the activations by least
% squares from the channels that remain instead: invert the REDUCED mixing
% matrix. That is the best linear estimate obtainable from a channel subset,
% and it is NOT the original decomposition -- hence the warning rather than
% a silent fallback.
    nComp = size(mixFull, 2);
    assert(sum(tf) > nComp, ...
        ['Only %d of %d ICA channels are present in the EEG, which is not more than ' ...
         'the %d components -- the activations are not identifiable. This usually ' ...
         'means the ICA and the EEG are not the same recording.'], ...
        sum(tf), numel(tf), nComp);
    warning('EEGCompareViewer:icaChanSubset', ...
        ['%d of %d ICA channels are absent from the EEG (%s).\n' ...
         'Re-projecting onto the %d channels present by least squares -- the ' ...
         'activations and the cleaned data are APPROXIMATE.'], ...
        sum(~tf), numel(tf), strjoin(labels(~tf), ', '), sum(tf));

    EEG.icawinv     = mixFull(tf, :);
    EEG.icasphere   = pinv(EEG.icawinv);   % least-squares unmixing over the retained channels
    EEG.icaweights  = eye(nComp);          % everything is folded into the sphere
    EEG.icachansind = idx(tf)';
end

function checkICAProvenance(file1, icafile)
% CHECKICAPROVENANCE An ICA is only meaningful for the recording it was
% trained on. A subject/session mismatch surfaces further down as "missing
% channels", which is a symptom rather than the problem -- so name it here.
    tok = @(s) unique(string(regexp(s, '(sub|ses)-[A-Za-z0-9]+', 'match')));
    if ~isequal(tok(file1), tok(icafile))
        warning('EEGCompareViewer:icaProvenance', ...
            ['ICA and EEG look like DIFFERENT recordings: EEG is [%s], ICA is [%s]. ' ...
             'Components from one recording do not describe another.'], ...
            strjoin(tok(file1), ' '), strjoin(tok(icafile), ' '));
    end
end

