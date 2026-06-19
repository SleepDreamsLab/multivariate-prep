function files = collectScoringFiles(scoringBase)
% JSON preferred, CSV fallback.
    f = dir(fullfile(scoringBase, '**', '*.json'));
    if isempty(f), f = dir(fullfile(scoringBase, '**', '*.csv')); end
    if isempty(f), files = {}; return; end
    files = fullfile({f.folder}, {f.name})';
end
