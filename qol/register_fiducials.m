function [chanlocs_reg, R, origin] = register_fiducials(chanlocs, varargin)
%register_fiducials Register channel positions to a fiducial-based head frame.
%
%   chanlocs_reg = register_fiducials(chanlocs) builds a head
%   coordinate frame from the Nasion, LPA and RPA fiducials found in an
%   EEGLAB chanlocs struct array (e.g. from readlocs() on a dome-solved
%   EGI .sfp file) and rotates/translates every channel into it:
%       +X = anterior  (toward nasion)
%       +Y = left      (toward LPA)
%       +Z = superior
%   This is the frame topoplot/headplot expect. Distorted/asymmetric
%   topographies caused by the raw photogrammetry frame not being
%   head-centered should resolve after this.
%
%   [chanlocs_reg, R, origin] = register_fiducials(chanlocs, ...)
%   also returns the 3x3 rotation matrix and 1x3 origin (both expressed
%   in the ORIGINAL coordinate frame), so you can apply the identical
%   registration to anything else computed in that native space (e.g. a
%   leadfield/head-model mesh):
%       xyz_new = (xyz_old - origin) * R';
%
%   Optional name-value arguments:
%     'NasionLabel'   - exact label of the nasion fiducial (default: auto-detect)
%     'LPALabel'      - exact label of the left fiducial   (default: auto-detect)
%     'RPALabel'      - exact label of the right fiducial  (default: auto-detect)
%     'KeepFiducials' - true/false, keep fiducial channels in the output
%                       (default: false, matching them being dropped before
%                       topoplot/ICA/etc.)
%
%   Auto-detection matches (case-insensitive) against:
%     Nasion: 'nasion','nas','nz','fidnz'
%     LPA:    'lpa','t9','fidt9'
%     RPA:    'rpa','t10','fidt10'
%
%   Example:
%     chanlocs     = readlocs('GSN256.sfp');
%     chanlocs_reg = register_fiducials(chanlocs);
%     topoplot(power_vector, chanlocs_reg);
%
%   If you need the same transform for a leadfield computed in the
%   native dome frame:
%     [chanlocs_reg, R, origin] = register_fiducials(chanlocs);
%     mesh_vertices_reg = (mesh_vertices - origin) * R';

p = inputParser;
addParameter(p, 'NasionLabel',   'FidNz',    @(x) ischar(x) || isstring(x));
addParameter(p, 'LPALabel',      'FidT9',    @(x) ischar(x) || isstring(x));
addParameter(p, 'RPALabel',      'FidT10',    @(x) ischar(x) || isstring(x));
addParameter(p, 'KeepFiducials', false, @islogical);
parse(p, varargin{:});
opt = p.Results;

labels = {chanlocs.labels};

nasion_idx = find_fiducial(labels, opt.NasionLabel, {'nasion','nas','nz','fidnz'}, 'Nasion');
lpa_idx    = find_fiducial(labels, opt.LPALabel,    {'lpa','t9','fidt9'},          'LPA');
rpa_idx    = find_fiducial(labels, opt.RPALabel,    {'rpa','t10','fidt10'},        'RPA');

nasion = get_xyz(chanlocs, nasion_idx, 'Nasion');
lpa    = get_xyz(chanlocs, lpa_idx,    'LPA');
rpa    = get_xyz(chanlocs, rpa_idx,    'RPA');

% --- build a right-handed head frame: X=anterior, Y=left, Z=superior ---
origin = (lpa + rpa) / 2;

x_hat = nasion - origin;
x_hat = x_hat / norm(x_hat);

y_raw = lpa - rpa;
y_raw = y_raw - dot(y_raw, x_hat) * x_hat;   % Gram-Schmidt: orthogonalize vs x_hat
y_hat = y_raw / norm(y_raw);

z_hat = cross(x_hat, y_hat);                  % unit length & orthogonal by construction

R = [x_hat; y_hat; z_hat];   % rows = new basis vectors, expressed in the original coords

% --- apply to every channel (fiducials included, in case they're kept) ---
chanlocs_reg = chanlocs;
for i = 1:numel(chanlocs)
    xyz_old = [chanlocs(i).X, chanlocs(i).Y, chanlocs(i).Z];
    if numel(xyz_old) ~= 3 || any(isnan(xyz_old))
        continue  % leave anything without valid coords untouched
    end
    xyz_new = (xyz_old - origin) * R';
    chanlocs_reg(i).X = xyz_new(1);
    chanlocs_reg(i).Y = xyz_new(2);
    chanlocs_reg(i).Z = xyz_new(3);
end

% --- drop fiducials unless requested otherwise ---
if ~opt.KeepFiducials
    chanlocs_reg([nasion_idx, lpa_idx, rpa_idx]) = [];
end

% --- recompute polar/spherical fields topoplot/headplot rely on ---
if exist('convertlocs', 'file')
    chanlocs_reg = convertlocs(chanlocs_reg, 'cart2all');
else
    warning('register_fiducials:noConvertlocs', ...
        ['EEGLAB function convertlocs.m not found on path. theta/radius/sph_*\n' ...
         'fields were NOT updated. Add EEGLAB to the path and run:\n' ...
         '  chanlocs = convertlocs(chanlocs, ''cart2all'');\n' ...
         'before plotting, or topoplot will use the stale angles.']);
end

% --- sanity-check printout ---
nas_new = (nasion - origin) * R';
lpa_new = (lpa    - origin) * R';
rpa_new = (rpa    - origin) * R';
fprintf('Registration check (expect Nasion~[+,~0,~0], LPA~[~0,+,~0], RPA~[~0,-,~0]):\n');
fprintf('  Nasion -> X=%6.2f  Y=%6.2f  Z=%6.2f\n', nas_new(1), nas_new(2), nas_new(3));
fprintf('  LPA    -> X=%6.2f  Y=%6.2f  Z=%6.2f\n', lpa_new(1), lpa_new(2), lpa_new(3));
fprintf('  RPA    -> X=%6.2f  Y=%6.2f  Z=%6.2f\n', rpa_new(1), rpa_new(2), rpa_new(3));

end

% ----------------------------------------------------------------------
function idx = find_fiducial(labels, userLabel, candidates, name)
if ~isempty(userLabel)
    idx = find(strcmpi(labels, userLabel));
    if numel(idx) ~= 1
        error('register_fiducials:badLabel', ...
            '%s: expected exactly one channel labeled ''%s'', found %d.', ...
            name, userLabel, numel(idx));
    end
    return
end
idx = find(ismember(lower(labels), candidates));
if numel(idx) ~= 1
    error('register_fiducials:autoDetectFailed', ...
        ['%s fiducial could not be auto-detected (found %d matches).\n' ...
         'Pass it explicitly, e.g.:\n' ...
         '  register_fiducials(chanlocs, ''%sLabel'', ''YourActualLabel'')'], ...
        name, numel(idx), name);
end
end

% ----------------------------------------------------------------------
function xyz = get_xyz(chanlocs, idx, name)
xyz = [chanlocs(idx).X, chanlocs(idx).Y, chanlocs(idx).Z];
if numel(xyz) ~= 3 || any(isnan(xyz))
    error('register_fiducials:missingCoords', ...
        '%s (channel ''%s'') has missing or NaN X/Y/Z coordinates.', ...
        name, chanlocs(idx).labels);
end
end