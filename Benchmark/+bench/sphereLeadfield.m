function G = sphereLeadfield(dipolePos, electrodePos, centre, radius)
% SPHERELEADFIELD  EEG leadfield of current dipoles inside a homogeneous conducting sphere.
%
%   G = bench.sphereLeadfield(dipolePos, electrodePos, centre, radius)
%
%   dipolePos     [nDip x 3] m, must lie inside the sphere
%   electrodePos  [nElec x 3] m, projected radially onto the sphere surface
%   centre        [1 x 3] m
%   radius        scalar m
%   G             [nElec x 3*nDip]; V = G * q for moments q = [qx1 qy1 qz1 qx2 ...].
%                 Conductivity and 1/(4*pi) are left out - callers scale to a target
%                 amplitude, so only the spatial pattern matters.
%
%   Closed-form single-shell solution (Zhang 1995, Phys Med Biol 40:335), the same
%   expression Brainstorm's bst_eeg_sph evaluates for one Berg dipole with mu = lambda = 1.
%   Used for the eyes, which are outside the cortical source space: a cortical leadfield
%   cannot represent them, and an infinite-medium dipole ignores the insulating scalp
%   boundary that shapes the frontal EOG field.

arguments
    dipolePos (:,3) double
    electrodePos (:,3) double
    centre (1,3) double
    radius (1,1) double {mustBePositive}
end

re = electrodePos - centre;
re = radius * re ./ vecnorm(re, 2, 2);
rq = dipolePos - centre;
if any(vecnorm(rq, 2, 2) >= radius)
    error('bench:sphereLeadfield:outside', 'Dipoles must lie inside the sphere.');
end

nElec = size(re, 1);
nDip  = size(rq, 1);
G = zeros(nElec, 3 * nDip);
for k = 1:nDip
    q = rq(k, :);
    qSq = sum(q.^2);
    dist = vecnorm(re - q, 2, 2);
    reDotQ = re * q';
    F  = dist .* (radius * dist + radius^2 - reDotQ);
    c1 = (2 * (reDotQ - qSq) ./ dist.^3 + 1 ./ dist - 1 / radius) / qSq;
    c2 = (2 ./ dist.^3 + (dist + radius) ./ (radius * F)) / qSq;
    G(:, 3*k-2:3*k) = (c1 - c2 .* reDotQ) * q + (c2 * qSq) .* re;
end
end
