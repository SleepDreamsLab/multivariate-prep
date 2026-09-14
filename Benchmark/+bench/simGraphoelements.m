function [S, events] = simGraphoelements(fwd, stageEpoch, srate, epochLength, classes, opts)
% SIMGRAPHOELEMENTS  Ground-truth sleep graphoelements projected through the cortical leadfield.
%
%   [S, events] = bench.simGraphoelements(fwd, stageEpoch, srate, epochLength, classes)
%   [S, events] = bench.simGraphoelements(..., 'templatebank', bank, 'epochorig', idx)
%
%   Inputs
%   ------
%   fwd           bench.buildForward output
%   stageEpoch    [1 x nEpochs] stage digit of each benchmark epoch
%   srate         Hz
%   epochLength   s
%   classes       struct array, see bench.defaultSignals
%
%   Name-value
%   ----------
%   epochorig     [1 x nEpochs] epoch index in the original recording (for the events table)
%   edgemargin    s kept free at both edges of every epoch - benchmark epochs are
%                 concatenated from non-adjacent parts of the night             (1.5)
%   guard         s between independent events (nested coupled events are exempt)  (0.3)
%   templatebank  struct with one field per class name, each .srate and .waveforms (cell of
%                 vectors). A class found there draws its time course from the bank instead
%                 of the parametric template - the hook for waveforms from a neural mass
%                 model. Placement, spatial projection and scaling are unchanged. (struct())
%   maxdelaybins  travelling waves are rendered with at most this many delay steps   (12)
%
%   Outputs
%   -------
%   S        [nChan x nSamples] single, uV, referenced to the recording reference
%   events   table, one row per event: id, trial_type, kind, onset (s), duration (s),
%            sample_start, sample_end, epoch, epoch_orig, stage, anchor, peakchan,
%            peakchan_idx, amplitude_uV, freq_Hz, speed_mps, coupled_to, template
%
%   Each event is one parametric (or bank) waveform emitted by a Gaussian-weighted patch of
%   normal-oriented cortical sources centred under an anchor channel. With a speed > 0 the
%   patch fires with delays along the anterior-posterior axis, so the field is not rank-1.
%   The event is scaled so that its average-referenced field reaches the drawn amplitude at
%   its peak channel; slow waves are sign-corrected to be surface-negative there.

arguments
    fwd struct
    stageEpoch (1,:) double
    srate (1,1) double {mustBePositive}
    epochLength (1,1) double {mustBePositive}
    classes struct
    opts.epochorig (1,:) double = []
    opts.edgemargin (1,1) double = 1.5
    opts.guard (1,1) double = 0.3
    opts.templatebank struct = struct()
    opts.maxdelaybins (1,1) double {mustBeInteger, mustBePositive} = 12
end

nChan        = size(fwd.Gctx, 1);
epochSamples = round(epochLength * srate);
nEpochs      = numel(stageEpoch);
nSamples     = nEpochs * epochSamples;
if isempty(opts.epochorig), opts.epochorig = 1:nEpochs; end

S        = zeros(nChan, nSamples, 'single');
occupied = false(1, nSamples);
margin   = round(opts.edgemargin * srate);
guard    = round(opts.guard * srate);
names    = {classes.name};
rows     = {};

%%% Long events first, while there is still room to place them
isSlow = strcmp({classes.kind}, 'slowwave');
for k = [find(isSlow), find(~isSlow)]
    c = classes(k);
    for iEpoch = 1:nEpochs
        iStage = find(c.stages == stageEpoch(iEpoch), 1);
        if isempty(iStage), continue; end
        rate = c.density(min(iStage, numel(c.density)));

        for iEvent = 1:poissonDraw(rate * epochLength / 60)
            [w, meta] = makeWaveform(c, srate, opts.templatebank);
            [E, info] = projectEvent(fwd, c, w, srate, opts.maxdelaybins);
            onset = tryPlace(occupied, iEpoch, epochSamples, size(E, 2), margin, guard);
            if isempty(onset), continue; end

            cols = onset:onset + size(E, 2) - 1;
            S(:, cols) = S(:, cols) + single(E);
            occupied(cols) = true;
            rows{end+1} = eventRow(numel(rows) + 1, c, cols, srate, iEpoch, ...
                opts.epochorig(iEpoch), stageEpoch(iEpoch), info, meta, 0); %#ok<AGROW>
            parentId = numel(rows);

            %%% Nested event, e.g. a spindle riding the up-state of a K-complex
            if isempty(c.coupleclass) || rand() >= c.coupleprob
                continue
            end
            nested = classes(strcmp(names, c.coupleclass));
            if isempty(nested)
                error('bench:simGraphoelements:unknownCoupleClass', ...
                    'Class %s nests unknown class %s.', c.name, c.coupleclass);
            end
            [w2, meta2] = makeWaveform(nested, srate, opts.templatebank);
            [E2, info2] = projectEvent(fwd, nested, w2, srate, opts.maxdelaybins);
            onset2 = onset + round((meta.negpeak + uniformDraw(c.couplelag)) * srate);
            if onset2 + size(E2, 2) - 1 > iEpoch * epochSamples - margin
                continue
            end
            cols2 = onset2:onset2 + size(E2, 2) - 1;
            S(:, cols2) = S(:, cols2) + single(E2);
            occupied(cols2) = true;
            rows{end+1} = eventRow(numel(rows) + 1, nested, cols2, srate, iEpoch, ...
                opts.epochorig(iEpoch), stageEpoch(iEpoch), info2, meta2, parentId); %#ok<AGROW>
        end
    end
end

if isempty(rows)
    events = struct2table(eventRow(0, classes(1), 1, srate, 1, 1, 0, emptyInfo(), ...
        struct('freq', NaN, 'negpeak', NaN, 'template', 0), 0), 'AsArray', true);
    events(1, :) = [];
else
    events = struct2table([rows{:}]', 'AsArray', true);
    events = sortrows(events, 'sample_start');
end
end

% -------------------------------------------------------------------------
function [w, meta] = makeWaveform(c, srate, bank)
% Unit waveform: spindles have an envelope peak of 1, slow waves a peak-to-peak of 1.
meta = struct('freq', NaN, 'negpeak', NaN, 'template', 0);

if isfield(bank, c.name) && ~isempty(bank.(c.name).waveforms)
    entry = bank.(c.name);
    meta.template = randi(numel(entry.waveforms));
    w = double(entry.waveforms{meta.template}(:)');
    if entry.srate ~= srate
        tOld = (0:numel(w) - 1) / entry.srate;
        w = interp1(tOld, w, 0:1/srate:tOld(end), 'pchip');
    end
    switch c.kind
        case 'spindle'
            w = w / max(abs(hilbert(w)));
        case 'slowwave'
            w = w / (max(w) - min(w));
            [~, iMin] = min(w);
            meta.negpeak = (iMin - 1) / srate;
        otherwise
            error('bench:simGraphoelements:kind', 'Unknown kind %s.', c.kind);
    end
    return
end

switch c.kind
    case 'spindle'
        duration = uniformDraw(c.duration);
        f0 = uniformDraw(c.freq);
        n  = max(2, round(duration * srate));
        u  = (0:n - 1) / (n - 1);
        %%% Hann envelope whose maximum sits at a random 35-60% of the event: warping u so
        %%% that u = peakAt maps to 0.5 keeps the peak value at exactly 1.
        peakAt = uniformDraw([0.35 0.6]);
        env = sin(pi * u.^(log(0.5) / log(peakAt))).^2;
        freq = f0 + uniformDraw(c.chirp) * (u - 0.5);
        w = env .* sin(2 * pi * cumsum(freq) / srate + 2 * pi * rand());
        meta.freq = f0;

    case 'slowwave'
        sdNeg = uniformDraw(c.negwidth);
        sdPos = uniformDraw(c.poswidth);
        tNeg  = 3.5 * sdNeg;
        tPos  = tNeg + uniformDraw(c.gap);
        t = (0:round((tPos + 3.5 * sdPos) * srate)) / srate;
        %%% Positive lobe amplitude sdNeg/sdPos makes both lobes enclose the same area, so
        %%% the event has no net DC - as a down-state followed by its up-state.
        w = -exp(-(t - tNeg).^2 / (2 * sdNeg^2)) + (sdNeg / sdPos) * exp(-(t - tPos).^2 / (2 * sdPos^2));
        w = w .* tukeywin(numel(t), 0.2)';
        w = w / (max(w) - min(w));
        meta.negpeak = tNeg;

    otherwise
        error('bench:simGraphoelements:kind', 'Unknown kind %s.', c.kind);
end
end

% -------------------------------------------------------------------------
function [E, info] = projectEvent(fwd, c, w, srate, maxBins)
info = emptyInfo();
info.anchor = c.anchors{randi(numel(c.anchors))};
src0 = bench.anchorSource(fwd, info.anchor);

dist  = vecnorm(fwd.GridLoc - fwd.GridLoc(src0, :), 2, 2);
patch = find(dist <= c.patchradius);
wts   = exp(-dist(patch).^2 / (2 * (c.patchradius / 2)^2));

info.speed = uniformDraw(c.speed);
if info.speed > 0
    %%% SCS +x points to the nasion, so -x is front-to-back propagation
    delays = (fwd.GridLoc(patch, :) - fwd.GridLoc(src0, :)) * [-1; 0; 0] / info.speed;
    shifts = round((delays - min(delays)) * srate);
    if numel(unique(shifts)) > maxBins
        step = max(shifts) / maxBins;
        shifts = round((min(floor(shifts / step), maxBins - 1) + 0.5) * step);
    end
else
    shifts = zeros(numel(patch), 1);
end
[steps, ~, group] = unique(shifts);

n = numel(w);
E = zeros(size(fwd.Gctx, 1), n + max(steps));
for b = 1:numel(steps)
    topo = fwd.Gctx(:, patch(group == b)) * wts(group == b);
    cols = steps(b) + 1:steps(b) + n;
    E(:, cols) = E(:, cols) + topo * w;
end

Eref = bench.avgRef(E);
switch c.kind
    case 'spindle'
        amp = max(abs(Eref), [], 2);
    case 'slowwave'
        amp = max(Eref, [], 2) - min(Eref, [], 2);
    otherwise
        error('bench:simGraphoelements:kind', 'Unknown kind %s.', c.kind);
end
[ampMax, info.peakchan_idx] = max(amp);
info.amplitude = uniformDraw(c.amplitude);
E = E * (info.amplitude / ampMax);

if strcmp(c.kind, 'slowwave') && sum(Eref(info.peakchan_idx, 1:n) .* w) < 0
    E = -E;   % surface-negative down-state at the peak channel
end
info.peakchan = fwd.labels{info.peakchan_idx};
end

% -------------------------------------------------------------------------
function info = emptyInfo()
info = struct('anchor', '', 'speed', 0, 'peakchan_idx', 1, 'peakchan', '', 'amplitude', NaN);
end

% -------------------------------------------------------------------------
function row = eventRow(id, c, cols, srate, iEpoch, epochOrig, stage, info, meta, parentId)
row = struct( ...
    'id',           id, ...
    'trial_type',   c.name, ...
    'kind',         c.kind, ...
    'onset',        (cols(1) - 1) / srate, ...
    'duration',     numel(cols) / srate, ...
    'sample_start', cols(1), ...
    'sample_end',   cols(end), ...
    'epoch',        iEpoch, ...
    'epoch_orig',   epochOrig, ...
    'stage',        stage, ...
    'anchor',       info.anchor, ...
    'peakchan',     info.peakchan, ...
    'peakchan_idx', info.peakchan_idx, ...
    'amplitude_uV', info.amplitude, ...
    'freq_Hz',      meta.freq, ...
    'speed_mps',    info.speed, ...
    'coupled_to',   parentId, ...
    'template',     meta.template);
end
