% TEST_BENCHMARK_UNITS  Unit checks for the +bench package on a toy spherical head.
%
%   Needs no recording and no GEDAI: a 64-channel montage and a radial-dipole leadfield are
%   written to tempdir, so bench.buildForward runs through its real file-reading path.
%   Covers the forward models, both simulators, both detectors, both scorers, background
%   synthesis and epoch selection. The GEDAI call itself (bench.runGedai) is exercised by
%   the integration run of bidsfun_gedai_benchmark.
%
%   Run from anywhere:  run('Benchmark/tests/test_benchmark_units.m')

set(0, 'DefaultFigureVisible', 'off');
testDir = fileparts(mfilename('fullpath'));
run(fullfile(testDir, '..', '..', 'dependancies.m'));
addpath(fullfile(testDir, '..'));
rng(7, 'twister');
srate = 250; epochLength = 30; epochSamples = epochLength * srate;
nPass = 0;

%% 1. Sphere leadfield: a dipole near the centre gives V = 3 q.r_hat / R^2 (times 4*pi*sigma)
R = 0.09;
elec = randn(50, 3); elec = R * elec ./ vecnorm(elec, 2, 2);
q = [0.3 -0.5 0.8];
G = bench.sphereLeadfield(1e-5 * [1 1 1], elec, [0 0 0], R);
expected = 3 * (elec / R) * q' / R^2;
assert(max(abs(G * q' - expected)) / max(abs(expected)) < 1e-3, 'sphere leadfield centre limit');
if exist('bst_eeg_sph', 'file') == 2
    dip = [0.02 -0.03 0.05];
    Gbst = bst_eeg_sph(dip, elec, [0 0 0], R, 1) * 4 * pi;
    Gours = bench.sphereLeadfield(dip, elec, [0 0 0], R);
    assert(max(abs(Gbst(:) - Gours(:))) / max(abs(Gours(:))) < 1e-4, 'sphere leadfield vs bst_eeg_sph');
end
nPass = nPass + 1; fprintf('PASS 1 sphere leadfield\n');

%% 2. Toy head files and bench.buildForward
[fwd, toy] = makeToyForward(tempdir);
assert(isequal(size(fwd.Gctx), [62, toy.nSrc]), 'Gctx size (two channels removed)');
assert(abs(fwd.sphere.radius - 0.09) < 1e-3, 'sphere fit radius');
src = bench.anchorSource(fwd, 'E1');
[~, nearest] = min(vecnorm(toy.sources - toy.elec(1, :) * 0.07 / 0.09, 2, 2));
assert(norm(toy.sources(src, :) - toy.sources(nearest, :)) < 0.02, 'anchor source lies under its channel');
eyeTopo = vecnorm(fwd.Keye, 2, 2);
[~, eyePeak] = max(eyeTopo);
assert(fwd.elecPos(eyePeak, 1) > 0.05, 'EOG field peaks at frontal electrodes');
nPass = nPass + 1; fprintf('PASS 2 buildForward / anchorSource\n');

%% 3. Graphoelements: amplitudes at the peak channel, margins, stage restriction
stageEpoch = [-2 -2 -2 -2 -3 -3 -3 -3 0 0];
classes = bench.defaultSignals();
for k = 1:numel(classes), classes(k).anchors = {'E1', 'E5', 'E9'}; end
[S, ev] = bench.simGraphoelements(fwd, stageEpoch, srate, epochLength, classes);
assert(height(ev) > 20, 'enough events injected');
assert(all(ismember(ev.stage, [-2 -3])), 'no events outside N2/N3');
Sref = bench.avgRef(S);
nChecked = 0;
for i = 1:height(ev)
    x = double(Sref(ev.peakchan_idx(i), ev.sample_start(i):ev.sample_end(i)));
    overlaps = ev.sample_start <= ev.sample_end(i) & ev.sample_end >= ev.sample_start(i);
    if nnz(overlaps) > 1, continue; end    % nested events add to each other's windows
    nChecked = nChecked + 1;
    switch ev.kind{i}
        case 'spindle'
            measured = max(abs(x));
        case 'slowwave'
            measured = max(x) - min(x);
    end
    assert(abs(measured / ev.amplitude_uV(i) - 1) < 0.02, 'event amplitude at peak channel');
end
assert(nChecked > 10, 'enough isolated events to check amplitudes');
epochStart = (ev.epoch - 1) * epochSamples + 1;
epochEnd = ev.epoch * epochSamples;
independent = ev.coupled_to == 0;
assert(all(ev.sample_start(independent) - epochStart(independent) >= 1.5 * srate), 'edge margin (start)');
assert(all(epochEnd - ev.sample_end >= 1.5 * srate - 1), 'edge margin (end)');
nPass = nPass + 1; fprintf('PASS 3 simGraphoelements (%d events)\n', height(ev));

%% 4. Artifacts: amplitudes and topography
cfg = bench.defaultArtifacts();
[A, evA] = bench.simArtifacts(fwd, stageEpoch, srate, epochLength, 'eog', cfg.eog);
Aref = bench.avgRef(A);
assert(height(evA) > 5, 'EOG events injected');
i = find(strcmp(evA.trial_type, 'eog_blink'), 1);
if ~isempty(i)
    x = Aref(evA.peakchan_idx(i), evA.sample_start(i):evA.sample_end(i));
    assert(abs(max(abs(x)) / evA.amplitude_uV(i) - 1) < 0.02, 'blink amplitude');
end
assert(mean(fwd.elecPos(unique(evA.peakchan_idx), 1)) > 0.03, 'EOG peak channels are frontal');

[A, evA] = bench.simArtifacts(fwd, stageEpoch, srate, epochLength, 'emg', cfg.emg);
Aref = bench.avgRef(A);
x = Aref(evA.peakchan_idx(1), evA.sample_start(1):evA.sample_end(1));
assert(abs(sqrt(mean(x.^2)) / evA.amplitude_uV(1) - 1) < 0.05, 'EMG RMS at peak channel');

[A, evA] = bench.simArtifacts(fwd, stageEpoch, srate, epochLength, 'ecg', cfg.ecg);
Aref = bench.avgRef(A);
assert(height(evA) == numel(stageEpoch), 'ECG in every epoch');
peakPerEpoch = max(abs(Aref(evA.peakchan_idx(1), :)));
assert(abs(peakPerEpoch / evA.amplitude_uV(1) - 1) < 0.15, 'ECG R-peak amplitude');
nPass = nPass + 1; fprintf('PASS 4 simArtifacts\n');

%% 5. Detectors: recall on clean injections, and frozen thresholds expose attenuation
nMin = 10; n = nMin * 60 * srate;
noise = pinkNoise(n, srate) * 10;
spindleTrue = zeros(1, n); swTrue = zeros(1, n);
onsets = round(linspace(5, nMin * 60 - 5, 40) * srate);
t = (0:srate - 1) / srate;
for o = onsets
    spindleTrue(o:o + srate - 1) = 30 * sin(pi * t).^2 .* sin(2 * pi * 13 * t);
    tw = (0:round(1.5 * srate)) / srate;
    wave = -exp(-(tw - 0.4).^2 / (2 * 0.12^2)) + 0.5 * exp(-(tw - 0.9).^2 / (2 * 0.24^2));
    swTrue(o + 2 * srate:o + 2 * srate + numel(tw) - 1) = 150 * wave / (max(wave) - min(wave));
end
x = noise + spindleTrue;
[dIn, thr] = bench.detectSpindles(x, srate);
recallIn = mean(arrayfun(@(o) any(dIn.onset <= o + srate & dIn.offset >= o), onsets));
dOut = bench.detectSpindles(noise + 0.3 * spindleTrue, srate, true(1, n), 'threshold', thr);
recallOut = mean(arrayfun(@(o) any(dOut.onset <= o + srate & dOut.offset >= o), onsets));
assert(recallIn > 0.9 && recallOut < 0.5, 'spindle detector recall %.2f -> %.2f', recallIn, recallOut);
x = noise + swTrue;
[dIn, thr] = bench.detectSlowWaves(x, srate);
troughs = onsets + 2 * srate + round(0.4 * srate);
tol = round(0.2 * srate);
recall = @(d) mean(arrayfun(@(tr) any(abs(d.negpeaksample - tr) <= tol), troughs));
recallIn = recall(dIn);
attenuated = noise + 0.1 * swTrue;
recallFrozen = recall(bench.detectSlowWaves(attenuated, srate, true(1, n), 'thresholds', thr));
recallRecalibrated = recall(bench.detectSlowWaves(attenuated, srate));
assert(recallIn > 0.9 && recallFrozen < 0.5 && recallFrozen <= recallRecalibrated, ...
    'slow-wave detector recall %.2f -> frozen %.2f, recalibrated %.2f', recallIn, recallFrozen, recallRecalibrated);
fprintf('  slow-wave recall %.2f -> %.2f frozen (%.2f recalibrated)\n', recallIn, recallFrozen, recallRecalibrated);
nPass = nPass + 1; fprintf('PASS 5 detectors\n');

%% 6. Signal scoring: perfect preservation and total removal
B = bench.avgRef(bench.simBackground(fwd, stageEpoch, srate, epochLength, 'synthetic', ...
    'stages', [-2 -3 0], 'targetrms', [20 30 15]));
rmsN2 = median(sqrt(mean(double(B(:, 1:4 * epochSamples)).^2, 2)));
assert(abs(rmsN2 / 20 - 1) < 0.15, 'synthetic background RMS %.1f', rmsN2);
m = bench.scoreSignals(Sref, B + Sref, B, B, ev, classes, stageEpoch, srate, epochLength);
v = @(tbl, cls, stg, name) tbl.value(strcmp(tbl.class, cls) & strcmp(tbl.stage, stg) & strcmp(tbl.metric, name));
% (B + S) - B in single precision leaves rounding error, so SER is large rather than Inf
assert(v(m, 'slowosc', 'all', 'SER_dB') > 60 && abs(v(m, 'slowosc', 'all', 'gain') - 1) < 1e-4, 'perfect preservation');
assert(v(m, 'slowosc', 'all', 'recall_in') == v(m, 'slowosc', 'all', 'recall_out'), 'identical recall');
m = bench.scoreSignals(Sref, B, B, B, ev, classes, stageEpoch, srate, epochLength);
assert(abs(v(m, 'spindle_fast', 'all', 'SER_dB')) < 1e-6 && abs(v(m, 'spindle_fast', 'all', 'gain')) < 1e-6, 'total removal');
assert(v(m, 'kcomplex', 'all', 'recall_out') <= v(m, 'kcomplex', 'all', 'recall_in'), 'removal cannot raise recall');
nPass = nPass + 1; fprintf('PASS 6 scoreSignals\n');

%% 7. Artifact scoring: perfect removal and no removal
[A, evA] = bench.simArtifacts(fwd, stageEpoch, srate, epochLength, 'eog', cfg.eog);
A = bench.avgRef(A);
stages = [-2 -3 0];
[PB, f] = bench.psdEpochs({B}, 1, stageEpoch, stages, srate, epochLength, 4);
psdB = struct('P', PB, 'f', f, 'stages', stages);
m = bench.scoreArtifacts(A, B, B, B, evA, 'eog', cfg.eog, stageEpoch, srate, epochLength, psdB, psdB);
assert(isinf(v(m, 'eog', 'all', 'ARR_dB')) && abs(v(m, 'eog', 'all', 'residual_gain')) < 1e-9, 'perfect removal');
assert(abs(v(m, 'eog', 'N2', 'psd_residual_excess')) < 1e-9, 'no residual excess');
m = bench.scoreArtifacts(A, B + A, B, B, evA, 'eog', cfg.eog, stageEpoch, srate, epochLength, psdB, psdB);
assert(abs(v(m, 'eog', 'all', 'ARR_dB')) < 1e-3 && abs(v(m, 'eog', 'all', 'residual_gain') - 1) < 1e-4, 'no removal');
assert(abs(v(m, 'eog', 'N2', 'psd_residual_excess') - 1) < 1e-3, 'full residual excess');
assert(v(m, 'eog', 'N2', 'psd_contamination_dB') > 0, 'EOG adds low-frequency power');
nPass = nPass + 1; fprintf('PASS 7 scoreArtifacts\n');

%% 8. Epoch selection and reading from an .fdt
nbchan = 8; nEp = 12;
raw = randn(nbchan, nEp * epochSamples, 'single');
raw(:, 3 * epochSamples + 1:4 * epochSamples) = raw(:, 3 * epochSamples + 1:4 * epochSamples) * 50; % epoch 4 is bad
fdt = fullfile(tempdir, 'bench_test.fdt');
fid = fopen(fdt, 'w'); fwrite(fid, raw, 'float32'); fclose(fid);
scoring = -2 * ones(1, nEp);
sel = bench.selectEpochs(fdt, nbchan, size(raw, 2), srate, epochLength, scoring, 'stages', -2, 'nepochs', 6, ...
    'candidatepool', 2, 'excludetransitions', false);
assert(~ismember(4, sel.epochs) && numel(sel.epochs) == 6 && issorted(sel.epochs), 'clean selection skips the bad epoch');
data = bench.readEpochs(fdt, nbchan, [2 5], epochSamples, 'demean', false);
assert(isequal(data, raw(:, [epochSamples + 1:2 * epochSamples, 4 * epochSamples + 1:5 * epochSamples])), 'readEpochs exact');
delete(fdt);
nPass = nPass + 1; fprintf('PASS 8 selectEpochs / readEpochs\n');

fprintf('\nAll %d benchmark unit checks passed.\n', nPass);

% =========================================================================
function [fwd, toy] = makeToyForward(folder)
% 64 electrodes on the upper part of a 9-cm sphere plus Cz, radial cortical dipoles at 7 cm
nElec = 64;
golden = pi * (3 - sqrt(5));
k = (0:nElec - 1)';
z = 1 - 1.3 * (k + 0.5) / nElec;               % down to ~30 degrees below the equator
r = sqrt(1 - z.^2);
elec = 0.09 * [r .* cos(golden * k), r .* sin(golden * k), z];
cz = [0 0 0.09];
nSrc = 800;
kk = (0:nSrc - 1)';
zs = 1 - 1.2 * (kk + 0.5) / nSrc;
rs = sqrt(1 - zs.^2);
sources = 0.07 * [rs .* cos(golden * kk), rs .* sin(golden * kk), zs];

Gain = bench.sphereLeadfield(sources, [elec; cz], [0 0 0], 0.09);
GridLoc = sources;
GridOrient = sources ./ vecnorm(sources, 2, 2);
lfFile = fullfile(folder, 'bench_toy_headmodel.mat');
save(lfFile, 'Gain', 'GridLoc', 'GridOrient');

sfpFile = fullfile(folder, 'bench_toy.sfp');
fid = fopen(sfpFile, 'w');
fprintf(fid, 'FidNz %.4f %.4f %.4f\n', 9.5, 0, 0);
fprintf(fid, 'FidT9 %.4f %.4f %.4f\n', 0, 7.5, 0);
fprintf(fid, 'FidT10 %.4f %.4f %.4f\n', 0, -7.5, 0);
for i = 1:nElec
    fprintf(fid, 'E%d %.4f %.4f %.4f\n', i, 100 * elec(i, :));
end
fprintf(fid, 'Cz %.4f %.4f %.4f\n', 100 * cz);
fclose(fid);

urlabels = arrayfun(@(i) sprintf('E%d', i), 1:nElec, 'uni', 0);
keep = setdiff(1:nElec, [10 20]);              % two channels "removed as bad"
chanlocs = struct('labels', urlabels(keep), 'urchan', num2cell(keep));
fwd = bench.buildForward(lfFile, sfpFile, chanlocs, urlabels, 'noteegchannels', nElec + 1:nElec + 20, ...
    'eyeoffset', [-2.5 3.2 -1.5], 'heartposition', [5 3 -32]);
toy = struct('elec', elec(keep, :), 'sources', sources, 'nSrc', nSrc);
toy.elec = elec;
delete(lfFile); delete(sfpFile);
end

% -------------------------------------------------------------------------
function x = pinkNoise(n, srate)
freqs = (0:n - 1) * srate / n;
freqs = min(freqs, srate - freqs);
shaping = 1 ./ max(freqs, 0.3);
shaping(freqs < 0.3) = 0;
x = real(ifft(fft(randn(1, n)) .* shaping));
x = x / std(x);
end
