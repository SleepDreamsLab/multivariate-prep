function lfCOV = loadrefcov(leadfielddir, p, nbchan, badChans, noteegchans)
    if isfolder(leadfielddir)
        lfFile = fullfile(leadfielddir, ['sub-' p.entities.sub], ['ses-' p.entities.ses], 'headmodel_surf_openmeeg.mat');
    else 
        lfFile = leadfielddir;
    end

    % Load leadfield matrix
    bstorm    = load(lfFile);

    % Remove not EEG chnnels
    noteegndx = intersect(1:size(bstorm.Gain,1), noteegchans);
    B         = bstorm.Gain;
    B(noteegndx, :) = [];

    % Remove bad channels
    if nbchan ~= size(B, 1)
        error('leadfield #chans is different from EEG #chans, potentially due to bad channel rejection indices')
    end
    % goodChans = setdiff(1:nbchan, find(removed_channels));
    % badChans  = intersect(1:nbchan, find(removed_channels));
    B(badChans, :) = [];

    % Compute leadfield covariance matrix
    B         = B - sum(B, 1) / (size(B, 1) + 1);
    lfCOV     = B * B';
end
