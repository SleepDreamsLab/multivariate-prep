clearvars; close all; clc

%%% preallocate
sEventsAll = {};
rEventsAll = {};

for ifile = 1:4

    %%% pathings
    eegfolder = '\\vs03\VS03-SandD-2\Sawtooth_detection\Data_collection\all';
    eegfile   = sprintf('H038_%d.mat', ifile);
    
    sawtwfolder = '\\vs03\VS03-SandD-2\Sawtooth_detection\Data_analysis\Data_processed\Presence_timestamps';
    sawtwfile = sprintf('H038_%d.csv', ifile);
    
    contrastfolder = '\\vs03\VS03-SandD-2\Sawtooth_detection\Data_analysis\Data_processed\Absence_timestamps';
    contrastfile = sprintf('H038_%d.csv', ifile);
    
    
    %%% read
    times1 = readtable(fullfile(sawtwfolder, sawtwfile));
    times2 = readtable(fullfile(contrastfolder, contrastfile));
    load(fullfile(eegfolder, eegfile));

    datavr = datavr - mean(datavr, 1);
    CZ = datavr(257, :);
    datavr(end, :) = [];
    
    sEvents = arrayfun(@(x) datavr(:, times1.start(x):times1.stop(x)), 1:size(times1,1), 'Uni', 0);
    rEvents = arrayfun(@(x) datavr(:, times2.start(x):times2.stop(x)), 1:size(times2,1), 'Uni', 0);

    sEventsAll = {sEventsAll{:}, sEvents{:}};
    rEventsAll = {rEventsAll{:}, rEvents{:}};
end

        opts.gedargs = {'contrast', 'spectral', 'peakfreq', 4.1, 'fwhm', 1.5, 'refmode', 'neighbour', 'neighbourdist', 1.5, 'segdur', 1, 'covnorm', 'trace', 'nperm', 2}
            G = ged(datavr,  'srate', 500, 'chanlocs',chanlocs,opts.gedargs{:});

    
G = ged(datavr, 'srate', 500, 'contrast', 'data', 'sdata', sEventsAll, 'rdata', rEventsAll, ...
        'covnorm', 'trace', 'nperm', 200, 'cvfolds', 3, 'plot', false, 'chanlocs',chanlocs);
plotged(G, 'ncomps', 4, 'xwindow', 30, 'acttype', 'signal')
plotged(G, 'ncomps', 4, 'xwindow', 60*5, 'acttype', 'wavelet_signal', 'cmap', slanCM('coolwarm', 20))

% figure; plot(datavr(257,:)); hold on
%%
comps = G.apply(datavr);
figure; 
plot(CZ(1, 1+500*30:(500*30+500*30))); hold on 
plot(comps(1,1+500*30:(500*30+500*30))); hold on 

% figure;
% plot(CZ(1, 1:500*30r) - comps(1,1:500*30 ))
