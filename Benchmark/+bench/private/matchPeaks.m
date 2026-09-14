function [hit, detIdx] = matchPeaks(evPeak, detPeak, tolerance)
% MATCHPEAKS  Match injected slow waves to detections by trough time.
%
%   hit(k)     true when a detected negative peak lies within tolerance samples of evPeak(k)
%   detIdx(k)  the nearest such detection (0 when none)
%
%   Stricter than window overlap: in N3 the background is dense with slow waves, and a
%   1-2 s window almost always overlaps some detection, whether or not the injected wave
%   survived.
nEvents = numel(evPeak);
hit = false(nEvents, 1);
detIdx = zeros(nEvents, 1);
if isempty(detPeak)
    return
end
detPeak = detPeak(:)';
for k = 1:nEvents
    [distance, j] = min(abs(detPeak - evPeak(k)));
    if distance <= tolerance
        hit(k) = true;
        detIdx(k) = j;
    end
end
end
