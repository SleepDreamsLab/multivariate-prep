function scoringFile = matchScoringFile(entities, allScoring)
% Match scoring file against all entities present in the EEG filename.
    scoringFile = '';
    if isempty(allScoring), return; end
    fields = fieldnames(entities);
    fields(contains(fields, 'acq')) = []; % ignore acq field
    mask   = true(size(allScoring));
    for k = 1:numel(fields)
        key  = fields{k};
        vals = regexp(allScoring, ['(?<=' key '-)[^_]+'], 'match', 'once', 'ignorecase');        
        mask = mask & contains(vals, entities.(key), 'IgnoreCase', true);
    end
    idx = find(mask, 1);
    if ~isempty(idx), scoringFile = allScoring{idx}; end
end
