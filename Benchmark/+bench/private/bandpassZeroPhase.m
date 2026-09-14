function y = bandpassZeroPhase(x, srate, band, order)
% BANDPASSZEROPHASE  Zero-phase Butterworth filter applied to each row of x.
%
%   band(1) <= 0 gives a low-pass, band(2) at or above ~Nyquist gives a high-pass. The filter
%   runs in second-order sections: the 0.5-Hz edge at 250 Hz is a small normalised frequency,
%   where transfer-function coefficients lose precision.
arguments
    x double
    srate (1,1) double
    band (1,2) double
    order (1,1) double = 4
end
nyquist = srate / 2;
if band(2) >= 0.98 * nyquist
    [z, p, k] = butter(order, band(1) / nyquist, 'high');
elseif band(1) <= 0
    [z, p, k] = butter(order, band(2) / nyquist, 'low');
else
    [z, p, k] = butter(order, band / nyquist, 'bandpass');
end
[sos, g] = zp2sos(z, p, k);
if isvector(x)
    y = reshape(filtfilt(sos, g, x(:)), size(x));
else
    y = filtfilt(sos, g, x')';
end
end
