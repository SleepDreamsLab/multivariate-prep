function name = buildSaveName(r, srate, modes, modesBB)
% BUILDSAVENAME  Build the run name describing a GEDAI configuration.
%
%   name = gedai.buildSaveName(r, srate, modes, modesBB)
%
%   modes/modesBB are the resolved per-stage-group GEDAIMode and GEDAIModeBB cell
%   arrays, so the name reflects what actually runs rather than a single nominal
%   value. Both are compacted by shortening 'auto' to 'a' and concatenated in
%   stage-group order; the broadband set is tagged 'BB'. For example
%   {'auto','auto','auto','auto+'} with {'auto-','auto-','auto-','auto+'} gives
%   '..._BBE10_aaaa+_BBa-a-a-a+_PT98_...'.

    name = sprintf( ...
        'ICA%s_BBonly%d_LCO%d_BBE%d_%s_BB%s_PT%d_BMT%d_ECT%d_SENSAI%d_KC%d_b1x%d_b2x%d_sr%d', ...
        r.ICAtype, r.broadbandOnly, r.GEDAILowCutOffFreq * 10, ...
        r.GEDAIBroadbandEpochSize, compactModes(modes), compactModes(modesBB), ...
        r.percentileThreshold, ...
        r.BBMinThreshold, r.GEDAIEnovaChannelThreshold, r.computeSENSAI, ...
        r.WeightKC * 100, round(r.boost1 * 10), round(r.boost2 * 10), srate);
end

function s = compactModes(modes)
% 'auto' -> 'a', 'auto+' -> 'a+', 'auto-' -> 'a-', concatenated in order.
    if ischar(modes) || isstring(modes)
        modes = cellstr(modes);
    end
    s = strjoin(cellfun(@(m) strrep(char(m), 'auto', 'a'), modes(:)', 'uni', 0), '');
end
