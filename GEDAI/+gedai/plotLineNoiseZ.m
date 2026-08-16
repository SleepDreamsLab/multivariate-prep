function plotLineNoiseZ(znoise, chanlocs, savefile, opts)
% PLOTLINENOISEZ  Topoplot of the residual line-noise z, the second detection round.
%
%   gedai.plotLineNoiseZ(znoise, chanlocs, savefile)
%   gedai.plotLineNoiseZ(znoise, chanlocs, savefile, threshold=4, title='sub-01', ...)
%
%   znoise     one value per entry of chanlocs, so pass the FULL montage with NaN for the
%              channels an earlier round already removed - they then show as absent rather
%              than as a hole the interpolation has to invent a value for, and the
%              electrode positions line up with the first-round topoplot.
%   chanlocs   EEG.urchanlocs, i.e. the montage as recorded
%   savefile   .png path
%
%   Companion to gedai.plotBadChannels, which draws the pre-Zapline round. Same statistic,
%   measured after line-noise removal instead of before, so the two are directly
%   comparable - and unlike the first round this one does remove channels, so everything
%   marked red here is missing from the saved data.
%
% See also: gedai.plotBadChannels, gedai.detectLineNoiseChannels

arguments
    znoise            double
    chanlocs          struct
    savefile          char
    opts.threshold    (1,1) double = 4
    opts.clim         (1,2) double = [0 10]
    opts.title        char         = ''
end

znoise = znoise(:);
nBad   = nnz(znoise > opts.threshold);

figure();
tiledlayout(1, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile();
drawTopoPanel(znoise, chanlocs, znoise > opts.threshold, opts.clim);
title({'Residual line noise, after Zapline/CleanLine', ...
       sprintf('%d channel(s) removed at z > %g', nBad, opts.threshold)});
if ~isempty(opts.title)
    subtitle(opts.title, 'Interpreter', 'none');
end

colormap('gray');
set(gcf, 'Color', 'w', 'Units', 'centimeters', 'Position', [2 2 11 11]);
gedai.printFigure(gcf, savefile);
close
end
