function EEG = assignChanlocs(EEG, BIDS, sfppath, sourceFile, p, fileID, opts)
% ASSIGNCHANLOCS  Attach channel locations to EEG from the SFP montage or BIDS ercp sidecars.
%   EEG = gedai.assignChanlocs(EEG, BIDS, sfppath, sourceFile, p, fileID)
%   EEG = gedai.assignChanlocs(..., 'forceoverwrite', true)
%
%   Skips reassignment when EEG.chanlocs already carries real coordinates - re-reading the
%   SFP would be wrong for a derivative that already had bad channels removed, since channel
%   k is then no longer the k-th SFP entry. Pass forceoverwrite to import anyway.
%
%   sfppath is tried first; the ercp sidecars ('_channels.tsv'/'_electrodes.tsv' next to
%   sourceFile) are the fallback, used only when sfppath is empty and BIDS.description.Name
%   is 'ercp'. p and fileID are the usual bids.internal.parse_filename output and the
%   entity-joined file identifier.
arguments
    EEG
    BIDS
    sfppath      char
    sourceFile   char
    p
    fileID       char
    opts.forceoverwrite (1,1) logical = false
end

hasChanlocs = isfield(EEG, 'chanlocs') && ~isempty(EEG.chanlocs) && ...
    isfield(EEG.chanlocs, 'X') && ~isempty([EEG.chanlocs.X]);

if hasChanlocs && ~opts.forceoverwrite
    fprintf('\nUsing channel locations already on the EEG struct\n')

elseif ~isempty(sfppath)
    sfpFile = gedai.matchSfpFile(sfppath, p.entities.sub, p.entities.ses);
    fprintf('\nReading %s ...\n', sfpFile)
    chanlocs     = readlocs(sfpFile);
    chanlocs_reg = register_fiducials(chanlocs);
    EEG.chanlocs = chanlocs_reg(1:EEG.nbchan);

    % Urchanlocs
    EEG.urchanlocs = EEG.chanlocs;
    for iCh = 1:numel(EEG.chanlocs)
        EEG.chanlocs(iCh).urchan = iCh;
    end

elseif strcmp(BIDS.description.Name, {'ercp'})
    chanfile = fullfile(fileparts(sourceFile), [fileID, '_channels.tsv']);
    elecfile = fullfile(fileparts(sourceFile), ['sub-' p.entities.sub, '_ses-' p.entities.ses, '_electrodes.tsv']);
    EEG = bids_importchanlocs(EEG, chanfile, elecfile);

    % Urchanlocs
    EEG.urchanlocs = EEG.chanlocs;
    for iCh = 1:numel(EEG.chanlocs)
        EEG.chanlocs(iCh).urchan = iCh;
    end

else
    % No coordinates available; caller proceeds without them.
end
end
