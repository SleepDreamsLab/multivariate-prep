function fwd = buildForward(lfFile, sfpFile, chanlocs, urlabels, opts)
% BUILDFORWARD  Forward models for signal and artifact injection.
%
%   fwd = bench.buildForward(lfFile, sfpFile, chanlocs, urlabels)
%   fwd = bench.buildForward(..., 'eyeoffset', [-2.5 3.2 -1.5], 'heartposition', [5 3 -32])
%
%   Everything is referenced to the recording reference (Cz): the potential at the reference
%   electrode is subtracted from every row. That is how the data on disk are referenced, and
%   bench.avgRef then applies the same average reference bidsfun_gedai does.
%
%   Inputs
%   ------
%   lfFile      Brainstorm head model (headmodel_surf_openmeeg.mat): Gain [nRows x 3*nSrc]
%               (or [nRows x nSrc] if already constrained), GridLoc, GridOrient. Rows are
%               the montage (E1..En) followed by the reference electrode, as used by
%               gedai.loadrefcov.
%   sfpFile     Electrode positions with fiducials (FidNz/FidT9/FidT10). Registered with
%               register_fiducials to +x anterior, +y left, +z up - the convention of the
%               Brainstorm SCS frame used for GridLoc.
%   chanlocs    EEG.chanlocs of the channels in the data (with .labels and .urchan)
%   urlabels    labels of EEG.urchanlocs, i.e. the full montage in leadfield row order
%
%   Name-value
%   ----------
%   noteegchannels  leadfield rows that are not montage channels            (257:300)
%   reflabel        label of the reference electrode in the SFP             ('Cz')
%   eyeoffset       cm, nasion to eye centre [anterior, lateral, superior]  ([-2.5 3.2 -1.5])
%   heartposition   cm, heart in the head frame                             ([5 3 -32])
%   eyefitradius    m, electrodes this close to the eyes define the eye sphere (0.10);
%                   the whole-head sphere is used when fewer than 30 qualify
%
%   Output fields
%   -------------
%   Gctx        [nChan x nSrc] normal-constrained cortical gain of the data channels
%   GcAll       [nUrChan x nSrc] the same for the full montage (anchor search)
%   peakRow, peakGain  per source, the montage row with the largest |gain| and its value
%   GridLoc     [nSrc x 3] m
%   labels, urlabels
%   elecPos     [nChan x 3] m, data channels; refPos [1 x 3] m
%   sphere      .centre, .radius (m), least-squares fit to all electrodes
%   eyeSphere   .centre, .radius (m), fit to the electrodes around the eyes
%   Keye        [nChan x 3] summed gain of both corneo-retinal dipoles (conjugate gaze),
%               homogeneous eye sphere (bench.sphereLeadfield)
%   Kheart      [nChan x 3] cardiac dipole gain, infinite homogeneous medium - the heart is
%               outside any head sphere, and at ~30 cm the scalp sees a smooth far field
%   eyePos, heartPos, nasion (m)

arguments
    lfFile char
    sfpFile char
    chanlocs struct
    urlabels cell
    opts.noteegchannels (1,:) double = 257:300
    opts.reflabel char = 'Cz'
    opts.eyeoffset (1,3) double = [-2.5 3.2 -1.5]
    opts.heartposition (1,3) double = [5 3 -32]
    opts.eyefitradius (1,1) double {mustBePositive} = 0.10
end

%%% Cortical leadfield
L = load(lfFile, 'Gain', 'GridLoc', 'GridOrient');
nRows = size(L.Gain, 1);
nSrc  = size(L.GridLoc, 1);
eegRows = setdiff(1:nRows, opts.noteegchannels);
if numel(eegRows) ~= numel(urlabels)
    error('bench:buildForward:rowMismatch', ...
        'Leadfield has %d montage rows but the recording montage has %d channels.', ...
        numel(eegRows), numel(urlabels));
end
refRow = setdiff(1:nRows, eegRows);

if size(L.Gain, 2) == 3 * nSrc
    %%% Constrain each source to its surface normal. Gain(:, 3j-2:3j) is source j's x/y/z.
    orient = L.GridOrient;
    Gc = L.Gain(:, 1:3:end) .* orient(:, 1)' + L.Gain(:, 2:3:end) .* orient(:, 2)' + ...
        L.Gain(:, 3:3:end) .* orient(:, 3)';
elseif size(L.Gain, 2) == nSrc
    Gc = L.Gain;
else
    error('bench:buildForward:gainSize', 'Gain has %d columns for %d sources.', size(L.Gain, 2), nSrc);
end
L = rmfield(L, 'Gain');

if isempty(refRow)
    %%% No reference row: the leadfield is taken to be reference-free already.
    GcAll = Gc(eegRows, :);
else
    GcAll = Gc(eegRows, :) - Gc(refRow(1), :);
end
clear Gc

goodRows = [chanlocs.urchan];
fwd.Gctx     = GcAll(goodRows, :);
fwd.GcAll    = GcAll;
[fwd.peakGain, fwd.peakRow] = max(abs(GcAll), [], 1);
fwd.GridLoc  = L.GridLoc;
fwd.labels   = {chanlocs.labels};
fwd.urlabels = urlabels(:)';

%%% Electrodes in the fiducial head frame
locs   = register_fiducials(readlocs(sfpFile), 'KeepFiducials', true);
labels = {locs.labels};
xyz    = [[locs.X]', [locs.Y]', [locs.Z]'];
xyz    = xyz * unitScale(xyz);

[found, idx] = ismember(fwd.labels, labels);
if ~all(found)
    error('bench:buildForward:missingElectrodes', 'Channels not in %s: %s', sfpFile, ...
        strjoin(fwd.labels(~found), ', '));
end
iRef = find(strcmpi(labels, opts.reflabel), 1);
iNas = find(ismember(lower(labels), {'fidnz', 'nasion', 'nas', 'nz'}), 1);
if isempty(iRef) || isempty(iNas)
    error('bench:buildForward:missingLandmark', ...
        'The SFP needs the reference electrode (%s) and a nasion fiducial.', opts.reflabel);
end
fwd.elecPos = xyz(idx, :);
fwd.refPos  = xyz(iRef, :);
fwd.nasion  = xyz(iNas, :);
if mean(fwd.elecPos(:, 3)) <= 0
    error('bench:buildForward:frame', ...
        'Electrodes sit below the fiducial plane on average - the SFP frame looks mirrored.');
end

%%% Head sphere through data channels and reference (EMG site directions)
allPos = [fwd.elecPos; fwd.refPos];
fwd.sphere = fitSphere(allPos);

%%% Eyes: corneo-retinal dipoles in a homogeneous sphere fitted locally. A whole-head sphere
%%% is centred well above the eyes (the net reaches far down the neck but not the face), and
%%% on DROP montages it puts the eye centres on or outside its surface. The electrodes
%%% within eyefitradius of the eyes describe the frontal curvature that shapes the EOG field.
offset = opts.eyeoffset / 100;
eyes = [fwd.nasion + [offset(1),  offset(2), offset(3)]; ...
        fwd.nasion + [offset(1), -offset(2), offset(3)]];
nearEyes = vecnorm(allPos - mean(eyes, 1), 2, 2) <= opts.eyefitradius;
if nnz(nearEyes) >= 30
    fwd.eyeSphere = fitSphere(allPos(nearEyes, :));
else
    fwd.eyeSphere = fwd.sphere;
end
eyeRadius = vecnorm(eyes - fwd.eyeSphere.centre, 2, 2);
maxRadius = 0.9 * fwd.eyeSphere.radius;
if any(eyeRadius > maxRadius)
    %%% The single-shell formula only holds inside the sphere. Pulling the eyes in along
    %%% the radius keeps their direction from the centre, which is what shapes the field.
    warning('bench:buildForward:eyesOutside', ...
        'Eye centres lie at %.2f of the eye-sphere radius; moved in to 0.9 along the radius.', ...
        max(eyeRadius) / fwd.eyeSphere.radius);
    eyes = fwd.eyeSphere.centre + (eyes - fwd.eyeSphere.centre) .* min(1, maxRadius ./ eyeRadius);
end
Keye = bench.sphereLeadfield(eyes, allPos, fwd.eyeSphere.centre, fwd.eyeSphere.radius);
Keye = Keye(1:end-1, :) - Keye(end, :);
fwd.Keye   = Keye(:, 1:3) + Keye(:, 4:6);
fwd.eyePos = eyes;

%%% Heart: current dipole in an infinite homogeneous medium
fwd.heartPos = opts.heartposition / 100;
d = allPos - fwd.heartPos;
Kheart = d ./ vecnorm(d, 2, 2).^3;
fwd.Kheart = Kheart(1:end-1, :) - Kheart(end, :);
end

% -------------------------------------------------------------------------
function sphere = fitSphere(pos)
% Algebraic least-squares sphere: |p|^2 = 2 p.c + (R^2 - |c|^2)
sol = [2 * pos, ones(size(pos, 1), 1)] \ sum(pos.^2, 2);
sphere.centre = sol(1:3)';
sphere.radius = sqrt(sol(4) + sum(sphere.centre.^2));
end

% -------------------------------------------------------------------------
function scale = unitScale(xyz)
% SFP files come in m, cm or mm depending on the exporter; a head radius of ~0.1 m tells.
headRadius = median(vecnorm(xyz - mean(xyz, 1), 2, 2));
if headRadius < 0.5
    scale = 1;
elseif headRadius < 50
    scale = 0.01;
else
    scale = 0.001;
end
end
