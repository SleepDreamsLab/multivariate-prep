function cmap = resolvecmap(spec, n)
% RESOLVECMAP  Turn a colormap option into an n-by-3 table.
%
%   Accepts the name of any colormap function ('hot', 'parula', 'turbo', ...)
%   or a ready-made n-by-3 matrix, so the plotters take a colormap the same way
%   MATLAB's own colormap() does. Shared by plotged and plotgednight.

arguments
    spec
    n (1,1) double {mustBePositive, mustBeInteger} = 256
end

if isnumeric(spec)
    if size(spec, 2) ~= 3
        error('resolvecmap:badMatrix', 'A colormap matrix must have three columns, not %d.', size(spec, 2));
    end
    cmap = spec;
    return
end

name = char(spec);
if exist(name, 'file') ~= 2 && exist(name, 'builtin') ~= 5
    error('resolvecmap:unknownColormap', ...
        '"%s" is not a colormap function on the path. Try hot, parula, turbo, gray, ...', name);
end
cmap = feval(name, n);
end
