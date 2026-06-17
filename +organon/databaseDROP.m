function Database = databaseDROP(BIDS, relativ_scoring_path, relative_save_path, relative_gps_path)
%CREATEBIDSDATABASE Create a table linking several files

arguments
    BIDS           % cell array of bids.layout objects
    relativ_scoring_path    = '../derivatives/scores/final'
    relative_save_path      = '../derivatives/exploratory-prep'
    relative_gps_path       = '../sourcedata/gps/**/solved'
end
%\sub-drop0001\solved\sub-drop0001_ses-t1_task-sleep_run-01_acq-domesolved_eeg.xml
% Pre-allocate database
VariableNames = {'FileID', 'EEG', 'Chanlocs', 'GPS', 'Scoring', 'Save'};
Database = table('Size', [0 numel(VariableNames)], ...
    'VariableNames', VariableNames, ...
    'VariableTypes', repelem({'cell'}, numel(VariableNames)));

% Loop through projects
for iBIDS = 1:numel(BIDS)
    BID = BIDS{iBIDS};    

    % Query sleep EEG files
    filesEEG = bids.query(BID, 'data', 'extension', '.vhdr', 'task', 'sleep');

    % Query channel/electrode files
    filesLOC = bids.query(BID, 'data', 'extension', '.tsv', 'task', 'sleep', ...
                          'suffix', {'channels', 'electrodes'});
    filesLOC = reshape(filesLOC, 2, [])';

    %%% --- Scoring files ---
    % Scoring files
    filesSCORING = dir(fullfile(BID.pth, relativ_scoring_path, '*eeg.csv'));
    filesSCORING = fullfile({filesSCORING.folder}, {filesSCORING.name})';

    %%% --- GPS files ---
    % GPS files
    filesGPS = dir(fullfile(BID.pth, relative_gps_path, '*eeg.sfp'));
    filesGPS = fullfile({filesGPS.folder}, {filesGPS.name})';
        

    % Match EEG files to scoring files
    [~, filesID, ~]     = fileparts(filesEEG);
    [~, scoringIDs, ~]  = fileparts(filesSCORING);
    [~, GPSIDs, ~]  = fileparts(filesGPS); 
    
    GPSIDs = strrep(GPSIDs, '_acq-domesolved', '');
    GPSIDs = strrep(GPSIDs, '-domesolved', '');
    
    [isMatch1, idxMatch1] = ismember(filesID, scoringIDs);
    [isMatch2, idxMatch2] = ismember(filesID, GPSIDs);
    isMatch = isMatch1 & isMatch2;


    % Only take those files where scoring exists
    filesID      = filesID(isMatch);
    filesEEG     = filesEEG(isMatch);
    filesLOC     = filesLOC(isMatch, :);
    filesSCORING = filesSCORING(idxMatch1(isMatch));
    filesGPS     = filesGPS(idxMatch2(isMatch));

    % Remove _eeg
    filesID = extractBefore(filesID, '_eeg');    

    % Where to save
    filesSAVE = char(java.io.File(fullfile(BID.pth, relative_save_path)).getCanonicalPath);
    filesSAVE = repelem({filesSAVE}, numel(filesID))';

    % Add rows to database
    newRows  = table(filesID, filesEEG, filesLOC, filesGPS, filesSCORING, filesSAVE, 'VariableNames', VariableNames);
    Database = [Database; newRows];
end
