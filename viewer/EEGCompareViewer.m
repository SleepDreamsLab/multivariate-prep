%% EEGCompareViewer
% Overlay two EEG recordings channel-by-channel for direct visual
% comparison, with epoch/stage navigation driven by a sleep scoring file.
% Requires EEGLAB (pop_interp, pop_select, ...) already on the path.
clc;
run('dependancies.m')

%% ===================== USER SETTINGS =====================
% file1       = '\\vs03.herseninstituut.knaw.nl\VS03-SandC-1\data\nin\data-drop\rawdata\derivatives\prep-ged\sub-drop0001\ses-t1\sub-drop0001_ses-t1_task-sleep_run-01_desc-zc2gedWakeBBAuto_eeg.set';   % first EEG file
% file1       = '\\vs03.herseninstituut.knaw.nl\VS03-SandC-1\data\nin\data-drop\rawdata\derivatives\prep-ged\sub-drop0001\ses-t1\sub-drop0001_ses-t1_task-sleep_run-01_desc-zc_eeg.vhdr';   % first EEG file
% file1       = '\\vs03.herseninstituut.knaw.nl\VS03-SandC-1\data\nin\data-drop\rawdata\derivatives\prep-ged\sub-drop0001\ses-t1\sub-drop0001_ses-t1_task-sleep_run-01_desc-zc2gedWakeBBAutoFSAutoPlus_eeg.set';   % second EEG file
% file2       = '\\vs03.herseninstituut.knaw.nl\VS03-SandC-1\data\nin\data-drop\rawdata\derivatives\prep-ged\sub-drop0001\ses-t1\sub-drop0001_ses-t1_task-sleep_run-01_desc-zc2gedWakeBBAutoPlusFSAutoPlus_eeg.set';   % second EEG file
file1       = '\\vs03.herseninstituut.knaw.nl\VS03-SandC-1\data\nin\data-drop\rawdata\derivatives\prep-ged\sub-drop0001\ses-t2\sub-drop0001_ses-t2_task-sleep_run-01_desc-zc2gedWakeBBAutoFSAutoPlus_eeg.set';   % second EEG file
scoringFile = '\\vs03.herseninstituut.knaw.nl\VS03-SandC-1\data\nin\data-drop\derivatives\scores\final\sub-drop0001_ses-t2_task-sleep_run-01_eeg.csv';                         % sleep scoring file, '' if none

net = 'EGI256'; % net type, passed to chans1020
targetrate = 100;
plotDecimation = 1; % only plot every Nth sample (display only; does not affect underlying data)

%% ===================== IMPORT + PREP =====================
EEG1 = fast_eeg_import(file1);

% ICA route
if targetrate < EEG1.srate; EEG1 = pop_resample(EEG1, targetrate); end
EEG = EEG1;

% root = '\\vs03.herseninstituut.knaw.nl\VS03-SandC-1\data\nin\data-drop\rawdata\derivatives\pamica\sub-drop0001_ses-t1_task-sleep_run-01_desc-zc2gedWakeBBAuto_eeg\zc-gedai-amica';
root = '\\vs03.herseninstituut.knaw.nl\VS03-SandC-1\data\nin\data-drop\rawdata\derivatives\pamica\sub-drop0001\ses-t1';
icafile = 'sub-drop0001_ses-t2_task-sleep_run-01_desc-zc2gedWakeBBAutoPlusFSAutoPlus2amica-Rej1Nmodel1F32Dll5e7MaxIter600Stride4_ica';
icafile = 'sub-drop0001_ses-t1_task-sleep_run-01_desc-noHP_ica'
pam = load(fullfile(root, [icafile '.mat']));   % for P (icasphere) + chanlabels
mod = loadmodout15(fullfile(root, icafile));   % W, S — byte-identical to Fortran

P = double(pam.icasphere);                        % k x nchan, the rank projection
W = double(mod.W(:,:,1));                         % explicit model index
S = double(mod.S(1:mod.num_pcs, 1:mod.data_dim)); % doc: S is num_pcs x data_dim

labels = string(pam.chanlabels(:));
assert(size(P,2) == numel(labels), 'P columns vs chanlabels mismatch');
assert(size(S,2) == size(P,1),     'amicaout rank vs P mismatch');

% amicaout lives in the P-space, so fold P into the sphere
mixFull   = pinv(W * S * P);                       % nchan_amica x ncomp, the mixing matrix
[tf, idx] = ismember(labels, string({EEG.chanlocs.labels}'));

% An ICA is only meaningful for the recording it was trained on. A
% subject/session mismatch surfaces further down as "missing channels",
% which is a symptom rather than the problem -- so name it here.
tok = @(s) unique(string(regexp(s, '(sub|ses)-[A-Za-z0-9]+', 'match')));
if ~isequal(tok(file1), tok(icafile))
    warning('EEGCompareViewer:icaProvenance', ...
        ['ICA and EEG look like DIFFERENT recordings: EEG is [%s], ICA is [%s]. ' ...
         'Components from one recording do not describe another.'], ...
        strjoin(tok(file1), ' '), strjoin(tok(icafile), ' '));
end

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
    % Some channels the decomposition was trained on are gone -- a different
    % derivative dropped them as bad. The exact unmixing W*S*P cannot be
    % formed, because it needs every column. Estimate the activations by
    % least squares from the channels that remain instead: invert the
    % REDUCED mixing matrix. That is the best linear estimate obtainable
    % from a channel subset, and it is NOT the original decomposition --
    % hence the warning rather than a silent fallback.
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

EEG = eeg_checkset(EEG, 'ica');

%% ===================== LABEL + SELECT =====================
EEG = iclabel(EEG);
EEG = selectcomps(EEG, 'ArtefactThreshold', .0, 'ManualQC', false, 'ICLabelClasses', [2:7]);
badComps = find(EEG.reject.gcompreject);
EEG2= pop_subcomp(EEG, badComps, 0);


scoringDigits = [];
if ~isempty(scoringFile)
    scoringDigits = scoreloader(scoringFile);
end

[~, fname1] = fileparts(file1);
file1 = fname1;
file2 = [fname1 ' + ICA'];

% [EEG1, chanmap] = chans1020(EEG1, false, 'add_eog', true, 'net', net, 'chanprefix', 'E');
% [EEG2, chanmap] = chans1020(EEG2, false, 'add_eog', true, 'net', net, 'chanprefix', 'E');
[EEG, chanmap] = chans1020(EEG, false, 'add_eog', true, 'net', net, 'chanprefix', 'E');

viewer.eegCompareViewer(EEG, [], scoringDigits, file1, file2, plotDecimation, fieldnames(chanmap));





