function [hit, detIdx] = matchEvents(evStart, evEnd, detOnset, detOffset, tolerance)
% MATCHEVENTS  Match injected events to detections by temporal overlap.
%
%   hit(k)     true when any detection overlaps [evStart(k) - tolerance, evEnd(k) + tolerance]
%   detIdx(k)  the detection with the largest overlap (0 when none)
%   All inputs in samples.
nEvents = numel(evStart);
hit = false(nEvents, 1);
detIdx = zeros(nEvents, 1);
if isempty(detOnset)
    return
end
detOnset = detOnset(:)';
detOffset = detOffset(:)';
for k = 1:nEvents
    overlap = min(evEnd(k) + tolerance, detOffset) - max(evStart(k) - tolerance, detOnset);
    [best, j] = max(overlap);
    if best >= 0
        hit(k) = true;
        detIdx(k) = j;
    end
end
end
