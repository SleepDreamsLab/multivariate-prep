function onset = tryPlace(occupied, iEpoch, epochSamples, len, margin, guard, maxTries)
% TRYPLACE  Random onset sample for an event of len samples inside epoch iEpoch.
%
%   The event keeps margin samples from both epoch edges - benchmark epochs are concatenated
%   from non-adjacent parts of the night, so an epoch edge is a discontinuity - and guard
%   samples from anything already marked in occupied. Returns [] when nothing fits.
arguments
    occupied (1,:) logical
    iEpoch (1,1) double
    epochSamples (1,1) double
    len (1,1) double
    margin (1,1) double
    guard (1,1) double
    maxTries (1,1) double = 50
end
onset = [];
epochStart = (iEpoch - 1) * epochSamples + 1;
epochEnd   = iEpoch * epochSamples;
lo = epochStart + margin;
hi = epochEnd - margin - len + 1;
if hi < lo
    return
end
for iTry = 1:maxTries
    candidate = randi([lo, hi]);
    a = max(epochStart, candidate - guard);
    b = min(epochEnd, candidate + len - 1 + guard);
    if ~any(occupied(a:b))
        onset = candidate;
        return
    end
end
end
