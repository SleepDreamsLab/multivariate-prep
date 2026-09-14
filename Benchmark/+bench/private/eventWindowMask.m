function mask = eventWindowMask(nSamples, startSample, endSample, pad)
% EVENTWINDOWMASK  Logical [1 x nSamples] mask covering every [start - pad, end + pad] window.
mask = false(1, nSamples);
for k = 1:numel(startSample)
    a = max(1, startSample(k) - pad);
    b = min(nSamples, endSample(k) + pad);
    mask(a:b) = true;
end
end
