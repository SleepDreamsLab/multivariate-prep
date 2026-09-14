function src = anchorSource(fwd, label)
% ANCHORSOURCE  Cortical source "under" a channel of the full montage.
%
%   src = bench.anchorSource(fwd, label)
%
%   Among the sources whose strongest (absolute, reference-subtracted) gain falls on that
%   channel, the one with the largest gain there. When no source peaks there, the source
%   with the highest gain at the channel relative to its own peak. Works from the leadfield
%   alone, so it needs no alignment between electrode and source coordinate frames, and it
%   also works for channels that were removed as bad (their leadfield rows remain).
arguments
    fwd struct
    label char
end
iRow = find(strcmp(fwd.urlabels, label), 1);
if isempty(iRow)
    error('bench:anchorSource:unknownChannel', 'Anchor channel %s is not in the montage.', label);
end
candidates = find(fwd.peakRow == iRow);
gainHere = abs(fwd.GcAll(iRow, :));
if isempty(candidates)
    [~, src] = max(gainHere ./ fwd.peakGain);
else
    [~, k] = max(gainHere(candidates));
    src = candidates(k);
end
end
