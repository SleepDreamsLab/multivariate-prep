function Database = databaseDIR(DIRECTORY, locFile, scoringext, relative_save_path)
%CREATEBIDSDATABASE Create a table linking several files

arguments
    DIRECTORY           % cell array of bids.layout objects
    locFile = [];
    scoringext = '_scoringfile.json'
    relative_save_path = 'TestFiles'
end

% Pre-allocate database
VariableNames = {'FileID', 'EEG', 'Chanlocs', 'Scoring', 'Save'};
Database = table('Size', [0 numel(VariableNames)], ...
    'VariableNames', VariableNames, ...
    'VariableTypes', repelem({'cell'}, numel(VariableNames)));

% Loop through projects
for iDIR = 1:numel(DIRECTORY)
    DIR = DIRECTORY(iDIR);    

    [~, filesID, ext]   = fileparts(DIR.name);
    filesEEG            = fullfile(DIR.folder, DIR.name);
    filesLOC            = locFile;
    % filesSCORING        = fullfile(DIR.folder, [filesID, '_scoringfile.json']);
    filesSCORING        = fullfile(DIR.folder, [filesID, scoringext]);
    filesSAVE           = char(java.io.File(fullfile(pwd, relative_save_path)).getCanonicalPath);

    % Add rows to database
    newRows  = table({filesID}, {filesEEG}, filesLOC, {filesSCORING}, {filesSAVE}, 'VariableNames', VariableNames);
    Database = [Database; newRows];
end
