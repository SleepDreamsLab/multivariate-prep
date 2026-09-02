function tf = hastopocoords(chanlocs)
% HASTOPOCOORDS  True when chanlocs carries coordinates topoplot can actually use.
%
%   Shared by plotged and plotgednight. A chanlocs that exists but holds only
%   labels - which happens with imported data whose montage was never resolved -
%   makes topoplot error out; without this check that one failure takes the whole
%   diagnostic figure with it, so both plotters degrade to a weight bar instead.

tf = false;
if isempty(chanlocs) || exist('topoplot', 'file') ~= 2
    return
end

%%% topoplot wants polar coordinates; it derives them from X/Y/Z when absent.
for f = {'theta', 'X'}
    if isfield(chanlocs, f{1})
        v = {chanlocs.(f{1})};
        if any(~cellfun(@isempty, v))
            tf = true;
            return
        end
    end
end
end
