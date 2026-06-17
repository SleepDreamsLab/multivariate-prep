function exploremodules(ModuleValues, Times_h, opts)
arguments
    ModuleValues struct
    Times_h 
    opts.Limits struct = struct()
    opts.Thresholds struct = struct()
    opts.nrow = 3;
    opts.ncol = 3;
    opts.colormap = 'copper';
end

% Default Limits
defaultLimits.SignalSTD                 = [0 5];
defaultLimits.VarianceSTD               = [0 .05];
defaultLimits.BigAmplitudes             = [0 1000];
defaultLimits.BigDifferences            = [0 600];
defaultLimits.BigDifferencesExtended    = [0 1000];
defaultLimits.MinCorrelation            = [0 1];
defaultLimits.MaxCorrelation            = [.998 1];
defaultLimits.SampleCorrelations        = [-1 1];
defaultLimits.SampleDifferences         = [0 20];
defaultLimits.SmoothMuscleSignalMax     = [0 800];

% Default Thresholds
defaultThresholds.SignalSTD             = 1.8;
defaultThresholds.VarianceSTD           = 0.01;
defaultThresholds.BigAmplitudes         = [-300 300];
defaultThresholds.BigDifferences        = 100;
defaultThresholds.BigDifferencesExtended= 300;
defaultThresholds.MinCorrelation        = .3;
defaultThresholds.MaxCorrelation        = .999;
defaultThresholds.SampleCorrelations    = .3;
defaultThresholds.SampleDifferences     = 5;
defaultThresholds.SmoothMuscleSignalMax = [0 0];

% Merge user-specified Limits/Thresholds with defaults
Limits = defaultLimits;
Thresholds = defaultThresholds;

userLimits = opts.Limits;
userThresholds = opts.Thresholds;

fn = fieldnames(userLimits);
for idot = 1:numel(fn)
    Limits.(fn{idot}) = userLimits.(fn{idot});
end

fn = fieldnames(userThresholds);
for idot = 1:numel(fn)
    Thresholds.(fn{idot}) = userThresholds.(fn{idot});
end

clf; 
set(gcf, 'color', 'w')

%%% Subfields
subfields = fieldnames(ModuleValues);

%%% Layout
layout = tiledlayout(opts.nrow, opts.ncol, ...
    'TileSpacing','compact', ...
    'Padding','compact');
cmap = colormap(opts.colormap);

%%% Colors for percentuals
prcmap = colormap('turbo');
prcmap = prcmap(1:floor(256/5):end, :);
prcmap = [
    0.0   0.45 0.70;   % minimum   - strong blue
    0.0   0.70 0.0;    % 25%       - strong green
    1.0   1.0  0.0;    % 50%       - yellow
    1.0   0.50 0.0;    % 75%       - orange
    0.85  0.0  0.0;    % maximum   - strong red
];


%%% Build figure
for ifield = 1:numel(subfields)
    
    nexttile(layout)
    hold on
    fname = subfields{ifield};

    data   = ModuleValues.(fname);
    nLines = size(data,1);

    h = plot(Times_h, data', '.', 'MarkerSize', 4);    
    arrayfun(@(i) set(h(i), 'Color', [cmap(i,:), .05]), 1:nLines);
    yline(Thresholds.(fname), 'r--')

    %%% Envelopes
    plot(Times_h, mean(data, 1), '.k', 'LineWidth', .1, 'MarkerSize', .5);    
    plot(Times_h, min(data, [], 1), '.b', 'LineWidth', .1, 'MarkerSize', .5, 'color', prcmap(1, :));
    plot(Times_h, prctile(data, 25, 1), '.c', 'LineWidth', .1, 'MarkerSize', .5, 'color', prcmap(2, :));     
    plot(Times_h, prctile(data, 50, 1), '.c', 'LineWidth', .1, 'MarkerSize', .5, 'color', prcmap(3, :));         
    plot(Times_h, prctile(data, 75, 1), '.y', 'LineWidth', .1, 'MarkerSize', .5, 'color', prcmap(4, :));        
    plot(Times_h, max(data, [], 1), '.r', 'LineWidth', .1, 'MarkerSize', .5, 'color', prcmap(5, :));

    xtickformat('%gh')
    title(fname)
    ylim(Limits.(fname))
    xlim([Times_h(1) Times_h(end)])
end

%%% Legend tyle
nexttile(layout)
hold on

% Define the labels for your percentiles
labels = {'Min', '25%', '50% (Median)', '75%', 'Max'};
    
% Plot dummy points for the Percentiles
for idot = 1:size(prcmap, 1)
    plot(nan, nan, '.', 'Color', prcmap(idot,:), 'MarkerSize', 15, 'DisplayName', labels{idot});
end
plot(nan, nan, '.k', 'MarkerSize', 10, 'DisplayName', 'Mean');
plot(nan, nan, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Threshold');

% Add the Legend to this specific tile
lgd = legend('Location', 'west', 'NumColumns', 1);
title(lgd, 'Legend')
lgd.Box = 'on';
axis off 

end
