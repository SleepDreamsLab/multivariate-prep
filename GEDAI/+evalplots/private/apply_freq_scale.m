function apply_freq_scale(ax, sc, logXLim, pow2ticks)
if strcmp(sc, 'log')
    set(ax, 'XScale', 'log', 'XLim', logXLim, 'XTick', pow2ticks);
end
end
