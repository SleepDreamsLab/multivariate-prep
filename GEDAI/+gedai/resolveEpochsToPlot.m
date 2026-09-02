function epochsToPlot = resolveEpochsToPlot(requested, scoringDigits)
    epochsToPlot = requested(:)';
    for score = -3:1
        idx = find(scoringDigits == score);
        if isempty(idx), continue; end
        %%% Five spread-out representatives of this stage. round(end/k) floors to 0 once
        %%% the stage has two epochs or fewer - a single epoch of one stage is common in a
        %%% night - so clamp to a valid index and let unique() below collapse the repeats.
        picks = [1, max(1, round(numel(idx) ./ [2 3 4 5]))];
        epochsToPlot = [epochsToPlot, reshape(idx(picks), 1, [])]; %#ok<AGROW>
        if score == -2 && numel(idx) >= 20
            epochsToPlot = [epochsToPlot, idx(5:5:20)]; %#ok<AGROW>
        end
    end
    epochsToPlot = unique(epochsToPlot);
end
