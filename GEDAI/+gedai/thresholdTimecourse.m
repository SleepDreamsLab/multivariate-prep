function T = thresholdTimecourse(EEG, opts)
% THRESHOLDTIMECOURSE  When and how much GEDAI's cleaning strength changed over a night.
%
%   T = gedai.thresholdTimecourse(EEG)
%   T = gedai.thresholdTimecourse(EEG, 'Scoring', scoringDigits, 'Plot', true)
%
% GEDAI records the artefact threshold it applied to every epoch of every band in
% EEG.etc.GEDAI.artifact_threshold_array_per_band, but each band has its own epoch
% length, so the arrays are not on a common time axis and are awkward to read. This puts
% them all on seconds and, optionally, plots them against the hypnogram.
%
% Why you would want this. With a sliding threshold (bidsfun_gedai's 'thresholdwindow')
% the cleaning strength is a continuous function of the data, which is the point - but it
% means the strength varies across a night, including across the stage transitions. If a
% transition is what you are analysing, a gradient in how much signal was removed is a
% confound you need to be able to see. This is how you see it: a flat threshold across
% your window of interest means the cleaning cannot be producing the effect; a ramp means
% you should check.
%
% Measured on one night, the broadband threshold rose from t = 3.3 about eight minutes
% before an awakening to t = 4.2 at the last sleep epoch. The cost in retained delta over
% that stretch was ~0.4 %, but the ramp is there and is worth inspecting per recording.
%
%   T   table with one row per band: name, epoch length, time vector (s), threshold
%       vector, and the Spearman correlation with sleep depth when Scoring is given.
%
% Name-value
% ----------
%   Scoring   vector of stage digits, one per 30-s epoch. Enables the hypnogram panel
%             and the correlation column. Use the ORIGINAL scoring, not the dilated
%             vector handed to GEDAI.
%   EpochLength  scoring epoch length in seconds (default 30).
%   Plot      draw the figure (default false).
%   Group     which entry of EEG.etc.GEDAI to read when a StageSpecific run produced
%             one per stage group (default 1). Whole-night runs have only one.
%
% See also bidsfun_gedai, gedai.wholeNightRuns.

arguments
    EEG
    opts.Scoring     double = []
    opts.EpochLength (1,1) double {mustBePositive} = 30
    opts.Plot        (1,1) logical = false
    opts.Group       (1,1) double {mustBeInteger, mustBePositive} = 1
end

if ~isfield(EEG, 'etc') || ~isfield(EEG.etc, 'GEDAI')
    error('gedai:thresholdTimecourse:noGEDAI', ...
        'EEG.etc.GEDAI is missing - was this file produced by GEDAI?');
end
G = EEG.etc.GEDAI;
if numel(G) < opts.Group
    error('gedai:thresholdTimecourse:noGroup', ...
        'EEG.etc.GEDAI has %d group(s); Group %d requested.', numel(G), opts.Group);
end
G = G(opts.Group);
if ~isfield(G, 'artifact_threshold_array_per_band')
    error('gedai:thresholdTimecourse:noThresholds', ...
        'This GEDAI output does not carry artifact_threshold_array_per_band.');
end

th = G.artifact_threshold_array_per_band;
nB = numel(th);

%%% Band epoch lengths: entry 1 is the broadband pass, the rest are the wavelet bands.
epochLen = zeros(1, nB);
epochLen(1) = G.broadband_epoch_size;
if nB > 1
    epochLen(2:nB) = G.epoch_sizes_per_wavelet_band(1:nB-1);
end

if isfield(G, 'freq_str_cell'), bandName = G.freq_str_cell(:);
else, bandName = arrayfun(@(b) sprintf('band %d', b), (1:nB)', 'uni', 0);
end

%%% Sleep depth as an ordinal, ordered the way artefact load is: N3 < N2 < REM < Wake.
depth = [];
if ~isempty(opts.Scoring)
    s = opts.Scoring(:)';
    depth = nan(size(s));
    depth(s == -3) = 1; depth(s == -2) = 2; depth(s == 0) = 3; depth(s == 1) = 4;
end

name = strings(nB,1); elen = zeros(nB,1);
time = cell(nB,1); thresh = cell(nB,1);
rho = nan(nB,1); tmin = nan(nB,1); tmax = nan(nB,1);
for b = 1:nB
    v = th{b}(:)';
    name(b)   = string(bandName{b});
    elen(b)   = epochLen(b);
    time{b}   = ((1:numel(v)) - 0.5) * epochLen(b);
    thresh{b} = v;
    tmin(b)   = min(v); tmax(b) = max(v);
    if ~isempty(depth) && numel(unique(v)) > 1
        ep = min(numel(depth), max(1, floor(time{b} / opts.EpochLength) + 1));
        d  = depth(ep);
        ok = ~isnan(d);
        if nnz(ok) > 2
            rho(b) = corr(v(ok)', d(ok)', 'type', 'Spearman');
        end
    end
end

T = table(name, elen, tmin, tmax, rho, time, thresh, ...
    'VariableNames', {'band','epochLength_s','t_min','t_max','rho_sleepDepth','time_s','threshold'});

if opts.Plot
    nPanel = min(nB, 8);
    figure('Color','w','Position',[60 60 1200 140*(nPanel+1)]);
    tl = tiledlayout(nPanel + double(~isempty(depth)), 1, ...
        'TileSpacing','compact','Padding','compact');
    tEnd = max(cellfun(@(x) x(end), time)) / 60;
    if ~isempty(depth)
        nexttile; stairs(((1:numel(depth))-0.5)*opts.EpochLength/60, depth, 'k', 'LineWidth', 1.1);
        set(gca,'YTick',1:4,'YTickLabel',{'N3','N2','REM','Wake'});
        ylim([0.5 4.5]); xlim([0 tEnd]); ylabel('stage'); box off
    end
    for b = 1:nPanel
        nexttile; plot(time{b}/60, thresh{b}, 'LineWidth', 1.2);
        xlim([0 tEnd]); grid on; box off; ylabel(name(b));
        if b == nPanel, xlabel('time (min)'); end
    end
    title(tl, 'GEDAI artefact threshold t applied over time, per band');
end
end
