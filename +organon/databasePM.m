function Database = databasePM(BIDS, opts)
%CREATEBIDSDATABASE Create a table linking several files

arguments
    BIDS           % cell array of bids.layout objects
    opts.task                    = {'Sleep', 'sleep'}    
    opts.recording               = '125Hz'
    opts.relative_scoring_path   = 'derivatives\scoring\scores\Manual_Checked'
    opts.relative_save_path      = 'derivatives\whatever'
end

% Pre-allocate database
VariableNames = {'FileID', 'EEG', 'Scoring', 'Chanlocs', 'Save'};
Database = table('Size', [0 numel(VariableNames)], ...
    'VariableNames', VariableNames, ...
    'VariableTypes', repelem({'cell'}, numel(VariableNames)));

% Loop through projects
for iBIDS = 1:numel(BIDS)
    BID = BIDS{iBIDS};    

    % Query sleep EEG files
    filesEEG = bids.query(BID, 'data', 'extension', '.vhdr', 'task', opts.task, 'recording', opts.recording);

    % % Query channel/electrode files
    % filesELEC = bids.query(BID, 'data', ...
    %     'extension', '.tsv', ...
    %     'suffix', {'electrodes'});
    % filesCHANS = bids.query(BID, 'data', ...
    %     'extension', '.tsv', ...
    %     'suffix', {'channels'});    
    % filesLOC = repmat({filesCHANS{1}, filesELEC{1}}, numel(filesEEG), 1);
    filesLOC = repmat({'', ''}, numel(filesEEG), 1);

    % Scoring files
    filesSCORING = dir(fullfile(BID.pth, opts.relative_scoring_path, '*.json'));
    filesSCORING = fullfile({filesSCORING.folder}, {filesSCORING.name})';
    % 
    % % Match EEG files to scoring files
    [~, filesID, ~]     = fileparts(filesEEG);
    [~, scoringIDs, ~]  = fileparts(filesSCORING);
    [isMatch, idxMatch] = ismember(filesID, scoringIDs);
    isMatch = 1:18; idxMatch = ones(1, 18, 'logical');
    
    % Only take those files where scoring exists
    filesID      = filesID(isMatch);
    filesEEG     = filesEEG(isMatch);
    filesLOC     = filesLOC(isMatch, :);
    filesSCORING = filesSCORING(idxMatch(isMatch));

    % Remove _eeg
    filesID = extractBefore(filesID, '_eeg');    

    % Where to save
    filesSAVE = char(java.io.File(fullfile(BID.pth, opts.relative_save_path)).getCanonicalPath);
    filesSAVE = repelem({filesSAVE}, numel(filesID))';

    % Add rows to database
    newRows  = table(filesID, filesEEG, filesSCORING, filesLOC, filesSAVE, 'VariableNames', VariableNames);
    Database = [Database; newRows];
end

end