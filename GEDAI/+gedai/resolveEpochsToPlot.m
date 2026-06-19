function epochsToPlot = resolveEpochsToPlot(requested, scoringDigits)
    epochsToPlot = requested(:)';
    for score = -3:1
        idx = find(scoringDigits == score);
        if isempty(idx), continue; end
        epochsToPlot = [epochsToPlot, idx(1), idx(round(end/2)), idx(round(end/3)), idx(round(end/4)), idx(round(end/5))]; %#ok<AGROW>
        if score == -2 && numel(idx) >= 20
            epochsToPlot = [epochsToPlot, idx(5:5:20)]; %#ok<AGROW>
        end
    end
    epochsToPlot = unique(epochsToPlot);
end
