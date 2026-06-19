function lfCOV = loadrefcov(leadfielddir, p, nbchan, removed_channels)
    if isfolder(leadfielddir)
        lfFile = fullfile(leadfielddir, ['sub-' p.entities.sub], ['ses-' p.entities.ses], 'headmodel_surf_openmeeg.mat');
    else 
        lfFile = leadfielddir;
    end
    bstorm    = load(lfFile);
    goodChans = setdiff(1:nbchan, find(removed_channels));
    B         = bstorm.Gain(goodChans, :);
    B         = B - sum(B, 1) / (size(B, 1) + 1);
    lfCOV     = B * B';
end
