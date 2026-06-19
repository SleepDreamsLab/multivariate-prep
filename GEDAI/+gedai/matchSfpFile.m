function sfp = matchSfpFile(CollectionRoot, SubjectName, SessName)
% Map a BIDS (subject, session) to its .sfp in the Data_collection tree.
%
% [VERIFY] Inferred from two (non-matching) examples:
%     sub-hpmam003  ->  H003_PM_AM           (trailing number fills the skeleton)
%     ses-S1        ->  S1
%     .../H###_PM_AM/SA_stim/S#/GPS/H###_PM_AM_GPS_S#_coordinates.sfp
% Edit this single function if your naming differs.
    lab    = erase(SubjectName, 'sub-');
    num    = regexp(lab, '\d+$', 'match', 'once');
    if isempty(num), sfp = ''; return; end
    subColl = sprintf('H%s_PM_AM', num);
    ses     = erase(SessName, 'ses-');
    gpsDir  = fullfile(CollectionRoot, subColl, 'SA_stim', ses, 'GPS');
    sfp     = fullfile(gpsDir, sprintf('%s_GPS_%s_coordinates.sfp', subColl, ses));
    if ~isfile(sfp)
        d = dir(fullfile(gpsDir, '*coordinates*.sfp'));
        if isempty(d), d = dir(fullfile(gpsDir, '*.sfp')); end
        if ~isempty(d), sfp = fullfile(d(1).folder, d(1).name); end
    end
end
