function [] = rasterplot(ArtefactsSamplewise, Srate, EpochLength, Scoring, opts)
arguments (Input)
    ArtefactsSamplewise
    Srate
    EpochLength
    Scoring = []
    opts.REM = 0;
    opts.Stagemap = dictionary( ...
        'W',    1, ...
        'REM',  0, ...
        'N1',   -1, ...
        'N2',   -2, ...
        'N3',   -3 ...
        );
end


%%% Downsample
nPnts               = size(ArtefactsSamplewise.Big, 2);
ArtefactsEpochWise  = sprep.resample_artefacts(ArtefactsSamplewise, Srate, EpochLength, nPnts);
fns                 = fieldnames(ArtefactsSamplewise);
nEpochs             = size(ArtefactsEpochWise.Big, 2);

%%% Times vector
EpochTimes = linspace(0, nEpochs*EpochLength/60/60, nEpochs);
PntsTimes  = linspace(0, nPnts/Srate/60/60, nPnts);

%%% Number of rows
nrow = numel(fns);
if ~isempty(Scoring)
    nrow = nrow + 1;
end

%%% Layout
ncol = 1;
layout = tiledlayout(nrow, ncol, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact')

%%% Scoring
if ~isempty(Scoring)
    nexttile()
    
    hold on
    stairs(EpochTimes, Scoring, 'k', 'LineWidth', 2)
    plot(EpochTimes(Scoring == opts.REM), opts.REM, 'r.')

    [SortedValues, Ndx] = sort(opts.Stagemap.values);
    SortedKeys = opts.Stagemap.keys;
    SortedKeys = SortedKeys(Ndx);
   
    yticks(SortedValues)
    yticklabels(SortedKeys)
    xlim([EpochTimes(1) EpochTimes(end)])
    xtickformat('%gh')
end

%%% Loop through submodules
for ifn = 1:numel(fns)
    fn = fns{ifn};

    submodule_data = ArtefactsEpochWise.(fn);
    nChans         = size(submodule_data, 1);    
    nPnts          = size(submodule_data, 2);

    nexttile()
    imagesc(PntsTimes, 1:nChans, ~submodule_data)
    title(fn)
    ylabel('Channel')
    xtickformat('%gh')
end
colormap("copper")

end