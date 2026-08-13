function topo_band_stage(PwrClean, PwrRaw, FreqsClean, FreqsRaw, StageScoring, Chanlocs, opts)
% GEDAI.EVALPLOTS.TOPO_BAND_STAGE  Band × sleep-stage topoplot grid, raw and clean.
%
%   Produces two figures (raw and clean), each a 6-band × 5-stage grid of
%   topoplots.  Rows: SWA, Theta, Alpha, Sigma, Beta, Gamma.
%   Columns: Wake, N1, N2, N3, REM.  Colour limits are per-band (derived
%   from clean data across all stages) and the same map as topo_band_power.

arguments
    PwrClean                    % chans × freqs × epochs
    PwrRaw
    FreqsClean   {mustBeVector}
    FreqsRaw     {mustBeVector}
    StageScoring {mustBeVector}
    Chanlocs
    opts.CLims    = []   % nBands×2 [lo hi] per band; NaN/missing rows → auto
    opts.SavePath = ''
end

bands = {
    'SWA 0.5–4 Hz',   [0.5  4];
    'Theta 4–8 Hz',   [4    8];
    'Alpha 8–12 Hz',  [8   12];
    'Sigma 12–16 Hz', [12  16];
    'Beta 16–30 Hz',  [16  30];
    'Gamma 30–45 Hz', [30  45];
};

stages = {
    'Wake',  1;
    'N1',   -1;
    'N2',   -2;
    'N3',   -3;
    'REM',   0;
};

stageOrder = [1, -1, -2, -3, 0];   % Wake, N1, N2, N3, REM
present    = unique(StageScoring);
present(present == 5) = 0;         % treat code 5 as REM
availableStages = stageOrder(ismember(stageOrder, present));

nBands  = size(bands,  1);
nStages = numel(availableStages);
nChans  = numel(Chanlocs);

% ── Compute mean log10 power: chans × bands × stages ─────────────────────
rawPwr   = nan(nChans, nBands, nStages);
cleanPwr = nan(nChans, nBands, nStages);

for iBand = 1:nBands
    fMaskRaw   = FreqsRaw   >= bands{iBand,2}(1) & FreqsRaw   <= bands{iBand,2}(2);
    fMaskClean = FreqsClean >= bands{iBand,2}(1) & FreqsClean <= bands{iBand,2}(2);
    if ~any(fMaskRaw) && ~any(fMaskClean), continue; end
    for iStage = 1:nStages
        sc    = availableStages(iStage);
        sMask = StageScoring == sc;
        if sc == 0, sMask = sMask | StageScoring == 5; end
        if ~any(sMask), continue; end
        if any(fMaskRaw)
            rawPwr(:,iBand,iStage)   = squeeze(mean(mean( ...
                log10(PwrRaw(:,fMaskRaw,sMask)   + eps), 2), 3));
        end
        if any(fMaskClean)
            cleanPwr(:,iBand,iStage) = squeeze(mean(mean( ...
                log10(PwrClean(:,fMaskClean,sMask) + eps), 2), 3));
        end
    end
end

% ── Per-band colour limits: user-supplied where provided, else data-driven ─
clims = nan(nBands, 2);
for iBand = 1:nBands
    if size(opts.CLims,1) >= iBand && all(isfinite(opts.CLims(iBand,:)))
        clims(iBand,:) = opts.CLims(iBand,:);
    else
        vals = cleanPwr(:,iBand,:);
        vals = vals(isfinite(vals));
        if ~isempty(vals)
            clims(iBand,:) = [floor(min(vals)*10)/10, ceil(max(vals)*10)/10];
        end
    end
end

warnState = warning('off', 'all');

% ── One figure for raw, one for clean ─────────────────────────────────────
datasets = {rawPwr,   'Raw EEG — band power per sleep stage',   'topo_band_stage_raw'; ...
            cleanPwr, 'Clean EEG — band power per sleep stage', 'topo_band_stage_clean'};

for iDs = 1:2
    pwr    = datasets{iDs,1};
    gtitle = datasets{iDs,2};
    suffix = datasets{iDs,3};

    fig = figure('Color','w','Units','centimeters', ...
        'Position',[2 2 nStages*6+3 nBands*5+1], 'Name', gtitle);
    tl  = tiledlayout(nBands, nStages, 'TileSpacing','compact','Padding','compact');
    title(tl, gtitle, 'FontSize',12,'FontWeight','bold');

    axGrid = gobjects(nBands, nStages);
    cbGrid = gobjects(nBands, 1);

    for iBand = 1:nBands
        cl = clims(iBand,:);
        for iStage = 1:nStages
            ax = nexttile;
            axGrid(iBand,iStage) = ax;
            data    = pwr(:,iBand,iStage);
            hasData = all(isfinite(data)) && all(isfinite(cl));

            if hasData
                topoplot(data, Chanlocs, 'maplimits', cl, ...
                    'electrodes','on','numcontour',0, 'conv', 'on');
                set(ax, 'CLim', cl);
            else
                axis(ax, 'off');
                text(0.5, 0.5, 'no data', 'Parent', ax, ...
                    'HorizontalAlignment','center','VerticalAlignment','middle', ...
                    'Units','normalized','FontSize',8,'Color',[0.6 0.6 0.6]);
            end

            % Column header on first band row
            if iBand == 1
                t = title(ax, stage_name(availableStages(iStage)), 'FontSize',10,'FontWeight','bold');
                t.Visible = 'on';
            end

            % Row label on first stage column — set Visible on so it shows
            % even when axis frame is hidden by topoplot
            if iStage == 1
                yl = ylabel(ax, bands{iBand,1}, 'FontSize',9,'FontWeight','bold');
                yl.Visible = 'on';
            end
        end

        % Colorbar attached to the last tile of each band row
        if all(isfinite(cl))
            cb = colorbar(axGrid(iBand,nStages));
            cb.FontSize = 7;
            cb.Limits   = cl;
            cb.Ticks    = linspace(cl(1), cl(2), 3);
            if iBand == nBands
                cb.Label.String = 'log_{10}(\muV^2/Hz)';
            end
            cbGrid(iBand) = cb;
        end
    end

    % Re-apply custom colormap per axis (topoplot resets figure-level map)
    for iBand = 1:nBands
        cl = clims(iBand,:);
        if any(~isfinite(cl)), continue; end
        for iStage = 1:nStages
            if isgraphics(axGrid(iBand,iStage))
                colormap(axGrid(iBand,iStage), custom_cmap());
                set(axGrid(iBand,iStage), 'CLim', cl);
            end
        end
        if isgraphics(cbGrid(iBand))
            cbGrid(iBand).Limits = cl;
            cbGrid(iBand).Ticks  = linspace(cl(1), cl(2), 3);
        end
    end

    set(fig, 'Color', 'w');
    save_fig(fig, opts.SavePath, suffix);
end

warning(warnState);
end
