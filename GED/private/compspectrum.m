function [pxx, hz] = compspectrum(x, srate)
% COMPSPECTRUM  Welch spectrum of a component time series.
%
%   Shared by plotged and plotgednight (hence private/). Falls back to a plain
%   Hann-tapered periodogram average when the Signal Processing Toolbox is
%   absent, so the diagnostic figures never depend on a toolbox.

x   = x(:);
win = min(numel(x), round(4 * srate));
if exist('pwelch', 'file') == 2
    [pxx, hz] = pwelch(x, hann(win), [], [], srate);
    return
end
nseg  = max(1, floor(numel(x) / win));
taper = 0.5 - 0.5 * cos(2 * pi * (0:win - 1)' / win);
pxx   = zeros(floor(win / 2) + 1, 1);
for s = 1:nseg
    seg = x((s - 1) * win + (1:win)) .* taper;
    amp = abs(fft(seg)).^2;
    pxx = pxx + amp(1:floor(win / 2) + 1);
end
pxx = pxx / (nseg * srate * sum(taper.^2));
hz  = linspace(0, srate / 2, numel(pxx))';
end
