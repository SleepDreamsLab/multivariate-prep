function sfp = matchSfpFile(bidsPath, subjectName, sessName)
% Dispatch to the study-specific SFP resolver based on fileID.
    if contains(subjectName, 'PM', 'IgnoreCase', true)
        sfp = sfpFromPM(bidsPath, subjectName, sessName);
    elseif contains(subjectName, 'DROP', 'IgnoreCase', true)
        sfp = sfpFromDrop(bidsPath, subjectName, sessName);
    else
        error('gedai:matchSfpFile:unknownStudy', ...
            'Cannot determine study from subject name "%s". Expected "PM" or "DROP".', subjectName);
    end
end

% -------------------------------------------------------------------------
function sfp = sfpFromPM(bidsPath, subjectName, sessName)
% PM: Data_collection is assumed to be a sibling of the BIDS folder. [VERIFY]
%     collRoot = fullfile(fileparts(bidsPath), 'Data_collection');
    collRoot = fullfile(bidsPath, 'sourcedata');
    lab = erase(subjectName, 'sub-');
    num = regexp(lab, '\d+$', 'match', 'once');
    if isempty(num), sfp = ''; return; end
    subColl = sprintf('H%s_PM_AM', num);   % [VERIFY skeleton]
    ses     = erase(sessName, 'ses-');
    gpsDir  = fullfile(collRoot, subColl, 'SA_stim', ses, 'GPS');
    sfp     = fullfile(gpsDir, sprintf('%s_GPS_%s_coordinates.sfp', subColl, ses));
    if ~isfile(sfp)
        d = dir(fullfile(gpsDir, '*coordinates*.sfp'));
        if isempty(d), d = dir(fullfile(gpsDir, '*.sfp')); end
        if ~isempty(d), sfp = fullfile(d(1).folder, d(1).name); end
    end
end

% -------------------------------------------------------------------------
function sfp = sfpFromDrop(bidsPath, subjectName, sessName)
% DROP: primary path is sourcedata/gps/<sub>/solved/*domesolved*.sfp;
% falls back to rawdata/<sub>/<ses>/eeg/*.sfp.
    solvedDir = fullfile(bidsPath, '..', 'sourcedata', 'gps', ['sub-' subjectName], 'solved');
    sessPat = ['*' sessName '*domesolved*.sfp'];
    d = dir(fullfile(solvedDir, sessPat));
    if isempty(d), d = dir(fullfile(solvedDir, '*.sfp')); end
    if ~isempty(d)
        sfp = fullfile(d(1).folder, d(1).name);
        return;
    end
end
