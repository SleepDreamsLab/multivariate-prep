function scoringFile = matchScoringFile(entities, allScoring)
% Match by sub + ses only (recording label may differ between EEG and scoring).
    scoringFile = '';
    if isempty(allScoring), return; end
    subs = regexp(allScoring, '(?<=sub-)[^_]+', 'match', 'once');
    sess = regexp(allScoring, '(?<=ses-)[^_]+', 'match', 'once');
    idx  = find(strcmp(subs, entities.sub) & strcmp(sess, entities.ses), 1);
    if ~isempty(idx), scoringFile = allScoring{idx}; end
end
