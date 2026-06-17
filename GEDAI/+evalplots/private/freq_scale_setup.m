function [scaleList, logXLim, pow2ticks] = freq_scale_setup(FreqLim, FreqScale)
switch FreqScale
    case 'linear', scaleList = {'linear'};
    case 'log',    scaleList = {'log'};
    case 'both',   scaleList = {'linear', 'log'};
    otherwise,     error('gedai.eval: FreqScale must be ''linear'', ''log'', or ''both''.');
end
logXLim   = [max(FreqLim(1), 1)  FreqLim(2)];
pow2ticks = [0, 2 .^ (0:10)];
pow2ticks = pow2ticks(pow2ticks >= logXLim(1) & pow2ticks <= logXLim(2));
end
