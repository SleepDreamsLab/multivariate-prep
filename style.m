% Colors ------------------------------------------------------------------
% *************************************************************************

% Dependancies
addpath('..\colormaps\Maps') % https://github.com/HuberSleepLab/colormaps

% Colormaps
colors = [];
colors.maps.div = colorcet('D1A');
colors.maps.lin = flipud(colorcet('L3')); % colorcet('L18');
colors.maps.cyc = colorcet('C3');
colors.maps.spec = spectral();
% colors.palette   = color_picker(10);
% colors.palette   = [colors.palette; colors.palette];

% Condition colors
colors.down = [0.0, 0.4, 0.95]; 
colors.sham = [0.2824 0.3294 0.3765];
colors.up   = [1.0000 0.2471 0.2039];

cyc = 50;
colors.stages = struct( ...
    'W', colors.maps.div(cyc*3,:), ...
    'N1', colors.maps.div(cyc*2,:), ...
    'N2', colors.maps.div(cyc*1,:), ...
    'N3', colors.maps.div(1,:), ...
    'REM', colors.maps.div(cyc*4,:), ...
    'None', [.1 .1 .1] ...
);