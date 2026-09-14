function [A, events] = simArtifacts(fwd, stageEpoch, srate, epochLength, type, cfg, opts)
% SIMARTIFACTS  Ground-truth EOG, EMG or ECG contamination with realistic topography.
%
%   [A, events] = bench.simArtifacts(fwd, stageEpoch, srate, epochLength, type, cfg)
%
%   Inputs
%   ------
%   fwd           bench.buildForward output (Keye, Kheart, elecPos, sphere)
%   stageEpoch    [1 x nEpochs] stage digit of each benchmark epoch
%   srate, epochLength
%   type          'eog', 'emg' or 'ecg'
%   cfg           the matching field of bench.defaultArtifacts
%
%   Name-value
%   ----------
%   epochorig     [1 x nEpochs] epoch index in the original recording
%   edgemargin    s kept free at the epoch edges (EOG, EMG)                      (1)
%   guard         s between events                                              (0.3)
%
%   Outputs
%   -------
%   A        [nChan x nSamples] single, uV, referenced to the recording reference
%   events   table with the columns of bench.simGraphoelements; EOG/EMG one row per event,
%            ECG one row per contaminated epoch (freq_Hz holds the heart rate in Hz)
%
%   Generators
%   ----------
%   EOG  Both eyes rotate together; the artifact is the change of each corneo-retinal
%        dipole, m(t) - m(0), with m = [cos v cos h, cos v sin h, sin v] for gaze angles
%        h (horizontal) and v (vertical), through fwd.Keye. Blinks are a vertical transient,
%        saccades a step-hold-return, slow eye movements a tapered sinusoidal drift.
%   EMG  Independent band-limited noise per channel plus a shared component, weighted by a
%        Gaussian spatial profile around a muscle site on the scalp. EMG reaches the
%        electrodes from muscle directly beneath them, so it is modelled in electrode space.
%   ECG  A far-field cardiac dipole (fwd.Kheart) driven by a three-lead vectorcardiogram,
%        continuous across all epochs of cfg.stages, heart rate with respiratory sinus
%        arrhythmia and beat-to-beat jitter.

arguments
    fwd struct
    stageEpoch (1,:) double
    srate (1,1) double {mustBePositive}
    epochLength (1,1) double {mustBePositive}
    type {mustBeMember(type, {'eog', 'emg', 'ecg'})}
    cfg struct
    opts.epochorig (1,:) double = []
    opts.edgemargin (1,1) double = 1
    opts.guard (1,1) double = 0.3
end

epochSamples = round(epochLength * srate);
nEpochs      = numel(stageEpoch);
if isempty(opts.epochorig), opts.epochorig = 1:nEpochs; end
A = zeros(size(fwd.Gctx, 1), nEpochs * epochSamples, 'single');

switch type
    case 'eog'
        [A, rows] = simEog(A, fwd, stageEpoch, srate, epochSamples, cfg, opts);
    case 'emg'
        [A, rows] = simEmg(A, fwd, stageEpoch, srate, epochSamples, cfg, opts);
    case 'ecg'
        [A, rows] = simEcg(A, fwd, stageEpoch, srate, epochSamples, cfg, opts);
    otherwise
        error('bench:simArtifacts:type', 'Unknown artifact type %s.', type);
end

if isempty(rows)
    events = struct2table(artifactRow(0, type, type, 1:2, srate, 1, 1, 0, '', 1, NaN, NaN), 'AsArray', true);
    events(1, :) = [];
else
    events = sortrows(struct2table([rows{:}]', 'AsArray', true), 'sample_start');
end
end

% =========================================================================
function [A, rows] = simEog(A, fwd, stageEpoch, srate, epochSamples, cfg, opts)
rows = {};
occupied = false(1, size(A, 2));
margin = round(opts.edgemargin * srate);
guard  = round(opts.guard * srate);
for subtype = {'blink', 'saccade', 'sem'}
    name = subtype{1};
    if ~isfield(cfg, name) || cfg.(name).density <= 0, continue; end
    sub = cfg.(name);
    for iEpoch = find(ismember(stageEpoch, cfg.stages))
        for iEvent = 1:poissonDraw(sub.density * epochSamples / srate / 60)
            [gaze, amplitude] = eyeCourse(name, sub, srate);
            E = fwd.Keye * gaze;
            onset = tryPlace(occupied, iEpoch, epochSamples, size(E, 2), margin, guard);
            if isempty(onset), continue; end
            [E, pc] = scaleToPeak(E, amplitude, 'maxabs');
            cols = onset:onset + size(E, 2) - 1;
            A(:, cols) = A(:, cols) + single(E);
            occupied(cols) = true;
            rows{end+1} = artifactRow(numel(rows) + 1, 'eog', ['eog_' name], cols, srate, ...
                iEpoch, opts.epochorig(iEpoch), stageEpoch(iEpoch), fwd.labels{pc}, pc, ...
                amplitude, NaN); %#ok<AGROW>
        end
    end
end
end

% -------------------------------------------------------------------------
function [gaze, amplitude] = eyeCourse(name, sub, srate)
% gaze: [3 x n] change of the unit corneo-retinal dipole from straight-ahead gaze
switch name
    case 'blink'
        duration = uniformDraw(sub.duration);
        n = round(duration * srate);
        t = (0:n - 1) / srate;
        tPeak = 0.35 * duration;
        bump = (t / tPeak).^3 .* exp(3 * (1 - t / tPeak));
        bump = bump .* tukeywin(n, 0.3)';
        v = deg2rad(10) * bump / max(bump);    % angle is irrelevant: scaled to amplitude
        h = zeros(1, n);
        amplitude = uniformDraw(sub.amplitude);

    case 'saccade'
        angle = uniformDraw(sub.angle);
        if rand() < sub.horizontalprob
            direction = pi * (rand() < 0.5);
        else
            direction = 2 * pi * rand();
        end
        nRamp = max(2, round(uniformDraw(sub.transition) * srate));
        ramp  = (1 - cos(pi * (0:nRamp - 1) / (nRamp - 1))) / 2;
        pad   = zeros(1, round(0.1 * srate));
        profile = [pad, ramp, ones(1, round(uniformDraw(sub.hold) * srate)), fliplr(ramp), pad];
        h = deg2rad(angle) * cos(direction) * profile;
        v = deg2rad(angle) * sin(direction) * profile;
        amplitude = angle * uniformDraw(sub.uvperdeg);

    case 'sem'
        n = round(uniformDraw(sub.duration) * srate);
        t = (0:n - 1) / srate;
        h = deg2rad(uniformDraw(sub.angle)) * sin(2 * pi * uniformDraw(sub.freq) * t + 2 * pi * rand()) ...
            .* tukeywin(n, 0.5)';
        v = zeros(1, n);
        amplitude = rad2deg(max(abs(h))) * uniformDraw(sub.uvperdeg);

    otherwise
        error('bench:simArtifacts:eogSubtype', 'Unknown EOG subtype %s.', name);
end
gaze = [cos(v) .* cos(h) - 1; cos(v) .* sin(h); sin(v)];
end

% =========================================================================
function [A, rows] = simEmg(A, fwd, stageEpoch, srate, epochSamples, cfg, opts)
rows = {};
occupied = false(1, size(A, 2));
margin = round(opts.edgemargin * srate);
guard  = round(opts.guard * srate);
nChan  = size(A, 1);

%%% Muscle sites: positions in the fiducial head frame (cm), each mapped to the nearest data
%%% channel. Fixed positions relative to nasion and ears stay put across nets; directions
%%% from a fitted sphere centre do not, because the centre moves with how far down the
%%% neck a net reaches.
siteChan = struct();
for s = cfg.sites
    if ~isfield(cfg.sitepositions, s{1})
        error('bench:simArtifacts:emgSite', 'EMG site %s has no entry in cfg.sitepositions.', s{1});
    end
    [~, siteChan.(s{1})] = min(vecnorm(fwd.elecPos - cfg.sitepositions.(s{1}) / 100, 2, 2));
end

band = [cfg.band(1), min(cfg.band(2), 0.45 * srate)];
for iEpoch = find(ismember(stageEpoch, cfg.stages))
    for iEvent = 1:poissonDraw(cfg.density * epochSamples / srate / 60)
        n = round(uniformDraw(cfg.duration) * srate);
        onset = tryPlace(occupied, iEpoch, epochSamples, n, margin, guard);
        if isempty(onset), continue; end

        site  = cfg.sites{randi(numel(cfg.sites))};
        sites = {site};
        if (startsWith(site, 'temporal') || startsWith(site, 'masseter')) && rand() < cfg.bilateralprob
            mirror = [site(1:end-1), char('L' + 'R' - site(end))];
            if isfield(siteChan, mirror), sites{end+1} = mirror; end %#ok<AGROW>
        end
        weight = zeros(nChan, 1);
        for s = sites
            d = vecnorm(fwd.elecPos - fwd.elecPos(siteChan.(s{1}), :), 2, 2);
            weight = weight + exp(-d.^2 / (2 * cfg.spread^2));
        end
        involved = find(weight > 0.02);

        %%% Filter with half a second of padding on both sides so the filter's edge
        %%% transients are cut away rather than injected
        pad = round(0.5 * srate);
        noise = bandpassZeroPhase(randn(numel(involved) + 1, n + 2 * pad), srate, band, 4);
        noise = noise(:, pad + 1:pad + n);
        noise = noise ./ sqrt(mean(noise.^2, 2));
        shared = uniformDraw(cfg.commonfraction);
        x = sqrt(1 - shared) * noise(1:end-1, :) + sqrt(shared) * noise(end, :);
        t = (0:n - 1) / srate;
        envelope = tukeywin(n, 0.4)' .* (1 + 0.3 * sin(2 * pi * uniformDraw([1 5]) * t + 2 * pi * rand()));

        E = zeros(nChan, n);
        E(involved, :) = weight(involved) .* x .* envelope;
        amplitude = uniformDraw(cfg.rms);
        [E, pc] = scaleToPeak(E, amplitude, 'rms');

        cols = onset:onset + n - 1;
        A(:, cols) = A(:, cols) + single(E);
        occupied(cols) = true;
        rows{end+1} = artifactRow(numel(rows) + 1, 'emg', ['emg_' strjoin(sites, '+')], cols, srate, ...
            iEpoch, opts.epochorig(iEpoch), stageEpoch(iEpoch), fwd.labels{pc}, pc, amplitude, NaN); %#ok<AGROW>
    end
end
end

% =========================================================================
function [A, rows] = simEcg(A, fwd, stageEpoch, srate, epochSamples, cfg, opts)
rows = {};
nSamples = size(A, 2);
t = (0:nSamples - 1) / srate;

%%% Wave model in cardiac phase (R peak at 0): P, Q, R, S, T. Row 1 is McSharry et al.
%%% (2003); rows 2-3 are illustrative orthogonal leads that make the loop rotate.
waveTheta = deg2rad([-70 -15 0 15 100]);
waveWidth = [0.25 0.1 0.1 0.1 0.4];
waveAmp   = [1.2 -5 30 -7.5 0.75; ...
             0.6  3 -12 -4  0.4; ...
             0.3 -2  6   2 -0.3];

%%% R-peak times: RR modulated by respiration (0.25 Hz) plus 2% beat-to-beat jitter
heartRate = uniformDraw(cfg.heartrate);
rr0 = 60 / heartRate;
rPeaks = -rand() * rr0;
while rPeaks(end) <= t(end) + rr0
    rPeaks(end+1) = rPeaks(end) + rr0 * (1 + cfg.rsa * sin(2 * pi * 0.25 * rPeaks(end)) + 0.02 * randn()); %#ok<AGROW>
end
beat  = discretize(t, rPeaks);
frac  = (t - rPeaks(beat)) ./ (rPeaks(beat + 1) - rPeaks(beat));
phase = 2 * pi * frac;
phase(frac > 0.5) = phase(frac > 0.5) - 2 * pi;

%%% Lead axes: the electrical axis, tilted at random, and two orthogonal leads
axis1 = cfg.axis(:)' / norm(cfg.axis);
tiltAxis = cross(axis1, randn(1, 3));
tiltAxis = tiltAxis / norm(tiltAxis);
axis1 = rotateVector(axis1, tiltAxis, deg2rad(uniformDraw([0, cfg.axisjitter])));
axis2 = cross(axis1, [0 0 1]);
if norm(axis2) < 1e-6, axis2 = cross(axis1, [1 0 0]); end
axis2 = axis2 / norm(axis2);
axis3 = cross(axis1, axis2);
leads = [axis1(:), axis2(:), axis3(:)];

phaseGrid = linspace(-pi, pi, 1000);
leadTemplate = vcg(phaseGrid, waveTheta, waveWidth, waveAmp);
baseline = mean(leadTemplate, 2);                 % no DC offset between stages

%%% Scale on one beat: the R-peak amplitude at the peak channel of the average-referenced field
amplitude = uniformDraw(cfg.amplitude);
[~, pc, scale] = scaleToPeak(fwd.Kheart * leads * (leadTemplate - baseline), amplitude, 'maxabs');

for iEpoch = find(ismember(stageEpoch, cfg.stages))
    cols = (iEpoch - 1) * epochSamples + 1:iEpoch * epochSamples;
    moment = leads * (vcg(phase(cols), waveTheta, waveWidth, waveAmp) - baseline);
    A(:, cols) = single(scale * fwd.Kheart * moment);
    rows{end+1} = artifactRow(numel(rows) + 1, 'ecg', 'ecg', cols, srate, iEpoch, ...
        opts.epochorig(iEpoch), stageEpoch(iEpoch), fwd.labels{pc}, pc, amplitude, heartRate / 60); %#ok<AGROW>
end
end

% -------------------------------------------------------------------------
function z = vcg(phase, waveTheta, waveWidth, waveAmp)
% [3 x n] lead voltages. Integrating McSharry's dz/dtheta = -a*dtheta*exp(-dtheta^2/(2b^2))
% over phase gives the closed form a*b^2*exp(-dtheta^2/(2b^2)) per wave.
z = zeros(size(waveAmp, 1), numel(phase));
for i = 1:numel(waveTheta)
    dTheta = mod(phase - waveTheta(i) + pi, 2 * pi) - pi;
    z = z + waveAmp(:, i) * (waveWidth(i)^2 * exp(-dTheta.^2 / (2 * waveWidth(i)^2)));
end
end

% -------------------------------------------------------------------------
function v = rotateVector(v, k, angle)
% Rodrigues rotation of v about unit axis k
v = v * cos(angle) + cross(k, v) * sin(angle) + k * dot(k, v) * (1 - cos(angle));
end

% =========================================================================
function [E, pc, scale] = scaleToPeak(E, amplitude, measure)
% Scale E so its average-referenced field reaches amplitude at the peak channel
Eref = bench.avgRef(E);
switch measure
    case 'maxabs'
        value = max(abs(Eref), [], 2);
    case 'rms'
        value = sqrt(mean(Eref.^2, 2));
    otherwise
        error('bench:simArtifacts:measure', 'Unknown measure %s.', measure);
end
[peak, pc] = max(value);
scale = amplitude / peak;
E = E * scale;
end

% -------------------------------------------------------------------------
function row = artifactRow(id, kind, name, cols, srate, iEpoch, epochOrig, stage, peakchan, pc, amplitude, freq)
row = struct( ...
    'id',           id, ...
    'trial_type',   name, ...
    'kind',         kind, ...
    'onset',        (cols(1) - 1) / srate, ...
    'duration',     numel(cols) / srate, ...
    'sample_start', cols(1), ...
    'sample_end',   cols(end), ...
    'epoch',        iEpoch, ...
    'epoch_orig',   epochOrig, ...
    'stage',        stage, ...
    'anchor',       '', ...
    'peakchan',     peakchan, ...
    'peakchan_idx', pc, ...
    'amplitude_uV', amplitude, ...
    'freq_Hz',      freq, ...
    'speed_mps',    0, ...
    'coupled_to',   0, ...
    'template',     0);
end
