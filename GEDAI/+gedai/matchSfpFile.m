function sfp = matchSfpFile(bidsPath, subjectName, sessName)
% MATCHSFPFILE  Resolve the SFP montage for one recording.
%   sfp = gedai.matchSfpFile(bidsPath, subjectName, sessName)
%
%   Dispatches to the study-specific resolver based on the subject name. Returns '' when
%   no montage exists for that subject/session - callers must treat that as a failure for
%   the recording, not as something to work around.
%
%   Every search below is scoped to the requested session. Electrode positions are
%   re-measured per session (the net is taken off and put back on between nights), so a
%   montage from another session of the same subject is wrong data, not an approximation:
%   it silently mislabels the geometry behind every topoplot and every interpolated
%   channel. There is deliberately no cross-session fallback.
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
% PM: sourcedata/<subColl>/SA_stim/<ses>/GPS/. The directory is already session-scoped,
%     so the widening searches inside it cannot cross into another session.
    sfp = '';
    collRoot = fullfile(bidsPath, 'sourcedata');
    lab = erase(subjectName, 'sub-');
    num = regexp(lab, '\d+$', 'match', 'once');
    if isempty(num), return; end
    subColl = sprintf('H%s_PM_AM', num);   % [VERIFY skeleton]
    ses     = erase(sessName, 'ses-');
    gpsDir  = fullfile(collRoot, subColl, 'SA_stim', ses, 'GPS');

    exact = fullfile(gpsDir, sprintf('%s_GPS_%s_coordinates.sfp', subColl, ses));
    if isfile(exact), sfp = exact; return; end

    d = dir(fullfile(gpsDir, '*coordinates*.sfp'));
    if isempty(d), d = dir(fullfile(gpsDir, '*.sfp')); end
    if isempty(d), return; end
    sfp = pickOne(d, subjectName, ses);
end

% -------------------------------------------------------------------------
function sfp = sfpFromDrop(bidsPath, subjectName, sessName)
% DROP: sourcedata/gps/<sub>/solved/*ses-<ses>*domesolved*.sfp, falling back to
%       rawdata/<sub>/ses-<ses>/eeg/*.sfp. The solved/ folder holds every session of the
%       subject side by side, so the session tag has to be in the pattern itself - this is
%       the search that used to widen to '*.sfp' and hand back a neighbouring session.
    sfp = '';
    ses       = erase(sessName, 'ses-');
    solvedDir = fullfile(bidsPath, '..', 'sourcedata', 'gps', ['sub-' subjectName], 'solved');

    d = dir(fullfile(solvedDir, ['*ses-' ses '_*domesolved*.sfp']));
    if isempty(d)
        %%% Same session, relaxed on the rest of the filename.
        d = dir(fullfile(solvedDir, ['*ses-' ses '_*.sfp']));
    end
    if isempty(d)
        rawDir = fullfile(bidsPath, ['sub-' subjectName], ['ses-' ses], 'eeg');
        d = dir(fullfile(rawDir, '*.sfp'));
    end
    if isempty(d), return; end
    sfp = pickOne(d, subjectName, ses);
end

% -------------------------------------------------------------------------
function sfp = pickOne(d, subjectName, ses)
% Take the first hit, but say so when the pattern was ambiguous. Several montages for one
% session means the naming has drifted, and which one you get is then down to sort order.
    sfp = fullfile(d(1).folder, d(1).name);
    if numel(d) > 1
        warning('gedai:matchSfpFile:ambiguous', ...
            '%d SFP files match sub-%s ses-%s; using "%s". Others: %s', ...
            numel(d), subjectName, ses, d(1).name, strjoin({d(2:end).name}, ', '));
    end
end
