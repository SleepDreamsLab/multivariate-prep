function value = uniformDraw(range)
% UNIFORMDRAW  One draw from U(range(1), range(2)). A scalar range is returned unchanged, so
% any simulation parameter can be given either as a fixed value or as a [min max] range.
if isscalar(range)
    value = range;
else
    value = range(1) + (range(2) - range(1)) * rand();
end
end
