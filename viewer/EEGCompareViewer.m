%% EEGCompareViewer
% Overlay two EEG recordings channel-by-channel for direct visual
% comparison, with epoch/stage navigation driven by a sleep scoring file.
% Requires EEGLAB (pop_interp, pop_select, ...) already on the path.
clc;

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
root = '\\vs03.herseninstituut.knaw.nl\VS03-SandC-1\data\nin\data-drop\rawdata\derivatives\pamica\sub-drop0001\ses-t2';
icafile = 'sub-drop0001_ses-t2_task-sleep_run-01_desc-zc2gedWakeBBAutoPlusFSAutoPlus2amica-Rej1Nmodel1F32Dll5e7MaxIter600Stride4_ica';
pam = load(fullfile(root, [icafile '.mat']));   % for P (icasphere) + chanlabels
mod = loadmodout15(fullfile(root, icafile));   % W, S — byte-identical to Fortran

P = double(pam.icasphere);                        % k x nchan, the rank projection
W = double(mod.W(:,:,1));                         % explicit model index
S = double(mod.S(1:mod.num_pcs, 1:mod.data_dim)); % doc: S is num_pcs x data_dim

labels = string(pam.chanlabels(:));
assert(size(P,2) == numel(labels), 'P columns vs chanlabels mismatch');
assert(size(S,2) == size(P,1),     'amicaout rank vs P mismatch');

% amicaout lives in the P-space, so fold P into the sphere
EEG.icaweights = W;
EEG.icasphere  = S * P;
EEG.icawinv    = pinv(EEG.icaweights * EEG.icasphere);   % same formula loadmodout15 uses

[tf, idx] = ismember(labels, string({EEG.chanlocs.labels}'));
assert(all(tf), 'channels missing from EEG: %s', strjoin(labels(~tf), ', '));
EEG.icachansind = idx(:)';

% cross-check the two export paths (scale-invariant, so normalization is moot)
C = abs(corr((EEG.icaweights*EEG.icasphere)', (double(pam.icaweights)*P)'));
r = max(C, [], 2);
fprintf('mod vs mat: median %.4f, min %.4f, n<0.99 = %d/%d\n', ...
        median(r), min(r), sum(r < 0.99), numel(r));

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





