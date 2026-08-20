function eegCompareViewer(EEG1, EEG2, scoringDigits, file1, file2, plotDecimation, chanLabelsToPlot)
% EEGCOMPAREVIEWER  Overlay one or two EEG recordings channel-by-channel for
% direct visual comparison, with epoch/stage navigation driven by a sleep
% scoring file. Requires EEGLAB (pop_interp, pop_select, ...) on the path.
%
%   eegCompareViewer(EEG1, EEG2, ...)   overlay two recordings
%   eegCompareViewer(EEG1, [],   ...)   open a single recording
%
% The single-recording form is meant for a recording that carries an ICA
% decomposition: the IC-subtracted version then becomes the second trace,
% so there is still something to compare against.
%
% All callbacks below are NESTED functions: they share this function's
% workspace directly (windowStart, diffMode, src, ...) instead of
% round-tripping state through get/setappdata on every redraw.
%
% plotDecimation: only every Nth sample is drawn (display-only). Default 5.
% chanLabelsToPlot: optional cell array of channel label strings to plot
% (matched against the channel set common to all inputs). Empty/omitted =
% plot all common channels. This is only the STARTING selection -- the
% "Channels..." button changes it at any time once the figure is up.
%
% ICA: if an input carries icaweights/icasphere, two extra controls light
% up (they stay greyed out otherwise):
%   IC view       -- plot component activations instead of channel traces
%                    (top 10 components by default).
%   Subtract ICs  -- overlay an EXTRA trace per ICA-bearing input with a
%                    chosen set of components projected out. Defaults to
%                    every non-brain component when ICLabel classifications
%                    are present. With two inputs that both carry an ICA
%                    decomposition this can put four traces on each row.
if nargin < 2, EEG2 = []; end
if nargin < 3, scoringDigits = []; end
if nargin < 4 || isempty(file1), file1 = 'signal 1'; end
if nargin < 5, file2 = ''; end
if nargin < 6 || isempty(plotDecimation), plotDecimation = 5; end
if nargin < 7, chanLabelsToPlot = {}; end

%%% --- Collect the inputs into a source list (1 or 2 recordings) ---
[~, n1] = fileparts(file1);
srcEEG  = {EEG1};
srcName = {n1};
if ~isempty(EEG2)
    if isempty(file2), file2 = 'signal 2'; end
    [~, n2]    = fileparts(file2);
    srcEEG{2}  = EEG2;
    srcName{2} = n2;
end
nSrc = numel(srcEEG);

%%% --- Align channels across the inputs (label match, stable order) ---
% The inputs can have had different bad channels dropped upstream, so their
% channel sets need not match -- keep only the intersection, and report
% exactly which labels were dropped (and from which input) so a mismatch
% reads as "these channels were bad in input X", not a silent reindex.
%
% NOTE: unlike earlier versions this does NOT reindex the data. Each
% source keeps its ORIGINAL channel rows, because icachansind indexes into
% those and would go stale under a reindex; the common-label set is mapped
% onto each source through src(k).rowOf instead.
labelsPerSrc = cellfun(@(E) {E.chanlocs.labels}, srcEEG, 'UniformOutput', false);
chLabelsFull = labelsPerSrc{1};
for k = 2:nSrc
    chLabelsFull = intersect(chLabelsFull, labelsPerSrc{k}, 'stable');
end
if isempty(chLabelsFull)
    error('eegCompareViewer:noCommonChans', 'The inputs share no channel labels.');
end
for k = 1:nSrc
    dropped = setdiff(labelsPerSrc{k}, chLabelsFull);
    if ~isempty(dropped)
        warning('eegCompareViewer:chanMismatch', ...
            'Using %d common channels; only in "%s" (dropped): %s.', ...
            numel(chLabelsFull), srcName{k}, strjoin(dropped, ', '));
    end
end
nChFull = numel(chLabelsFull);

srate = srcEEG{1}.srate;
nSamp = min(cellfun(@(E) size(E.data, 2), srcEEG));
for k = 2:nSrc
    if srcEEG{k}.srate ~= srate
        warning('eegCompareViewer:srateMismatch', ...
            'Sampling rates differ (%.3f vs %.3f Hz); using "%s"''s rate.', ...
            srate, srcEEG{k}.srate, srcName{1});
    end
end

%%% --- Build the per-source record ---
src = struct('name', {}, 'data', {}, 'rowOf', {}, 'hasICA', {}, 'unmix', {}, ...
             'winv', {}, 'icachansind', {}, 'posInICA', {}, 'nComp', {}, ...
             'icTick', {}, 'icText', {}, 'icClass', {}, 'classNames', {}, ...
             'subIdx', {}, 'subP', {});
for k = 1:nSrc
    E = srcEEG{k};
    src(k).name  = srcName{k};
    src(k).data  = single(E.data(:, 1:nSamp));   % original channel order, see note above
    [~, rowOf]   = ismember(chLabelsFull, labelsPerSrc{k});
    src(k).rowOf = rowOf;
    src(k)       = attachICA(src(k), E);
end
icaSrc = find([src.hasICA]);
if isempty(icaSrc)
    icaPrimary = [];
else
    icaPrimary = icaSrc(1);
end

%%% --- Which channels are plotted at startup ---
% A requested label can be legitimately absent -- the inputs may have had
% different channels dropped as bad upstream, so a label missing from the
% common set doesn't mean the caller made a mistake. Skip those with a
% warning; only error if NONE of the requested labels are available.
if isempty(chanLabelsToPlot)
    initChanIdx = 1:nChFull;
else
    [found, chanIdx] = ismember(chanLabelsToPlot, chLabelsFull);
    if ~any(found)
        error('eegCompareViewer:badChanLabel', ...
            'None of the requested channel labels are present in the common channel set: %s', ...
            strjoin(chanLabelsToPlot, ', '));
    end
    if ~all(found)
        warning('eegCompareViewer:chanLabelMissing', ...
            'Requested channel label(s) not in the common channel set (likely a bad channel dropped in one of the inputs) and will be skipped: %s', ...
            strjoin(chanLabelsToPlot(~found), ', '));
    end
    initChanIdx = sort(chanIdx(found));   % ascending -- keeps chLabelsFull's stable order
end

totalDur      = nSamp / srate;         % seconds
scoreEpochSec = 30;                    % fixed: sleep-scoring epoch length, used for epoch/stage alignment
nEpochs       = max(1, floor(totalDur / scoreEpochSec));

%%% --- Colours ---
% A source and its IC-subtracted overlay land on the SAME row, one on top
% of the other, so that pair has to carry the strongest contrast in the
% palette -- each raw colour is therefore paired with its rough complement
% (blue/orange, crimson/green) rather than with a neighbouring hue. A
% desaturated partner was tried first and was unreadable: at 0.5 pt the two
% lines were near-indistinguishable.
royalBlue    = [ 65 105 225] / 255;    % source 1, raw
orangeColor  = [235 125  10] / 255;    % source 1, ICs removed
redColor     = [220  20  60] / 255;    % source 2, raw -- crimson
greenColor   = [ 15 145  70] / 255;    % source 2, ICs removed
diffColor    = [120  45 190] / 255;    % violet -- a derived trace, distinct from all four
zeroLineGray = [0.75 0.75 0.75];
colRaw   = {royalBlue,   redColor};
colClean = {orangeColor, greenColor};
maxSeries = 2 * nSrc;                  % raw + cleaned per source

%%% --- Row spacing for the stacked traces (robust, subsampled) ---
% Deliberately computed from the FULL channel / component set, not just the
% currently-plotted rows, so the vertical scale doesn't jump around every
% time the selection changes. Channel view and IC view get their own
% spacing: component activations are in arbitrary units, typically orders
% of magnitude away from microvolts.
probeIdx    = 1:max(50, floor(nSamp / 20000)):nSamp;
probe       = cell(1, nSrc);
for k = 1:nSrc
    probe{k} = src(k).data(src(k).rowOf, probeIdx);
end
spacingChan = robustSpacing(probe);

spacingIC = spacingChan;
if ~isempty(icaSrc)
    probeAct = cell(1, numel(icaSrc));
    for j = 1:numel(icaSrc)
        k = icaSrc(j);
        probeAct{j} = src(k).unmix * src(k).data(src(k).icachansind, probeIdx);
    end
    spacingIC = robustSpacing(probeAct);
end

%%% --- Sleep-stage lookup (qol/scoreloader.m convention) ---
stageNames  = {'N3', 'N2', 'N1', 'REM', 'Wake'};
stageValues = [-3, -2, -1, 0, 1];

%%% --- Navigation / display state (shared with all nested callbacks) ---
windowStart   = 0;
displaySec    = scoreEpochSec;   % how much time is shown at once; user-adjustable, defaults to 30 s
diffMode      = false;
cleanOnly     = false;    % show only the last (cleanest) trace; mutually exclusive with diffMode
subtractOn    = false;    % overlay the IC-subtracted trace(s)
viewMode      = 'chan';   % 'chan' = channel traces, 'ic' = component activations
ampScale      = 1;        % vertical gain applied to the traces before offsetting
lastStateKey  = '';       % guards the legend/button refresh -- see redraw()
gridStepCandidates = [1 2 5 10 15 30 60 120 300 600 900 1800 3600 7200];   % "nice" seconds for the vertical time grid
maxGridLines       = 60;  % cap on vertical gridlines so a big window doesn't paint hundreds of them

%%% --- Row-selection state ---
% Everything below is (re)built by applyRowSelection() -- once at startup,
% and again whenever the view mode or the channel/IC selection changes.
% Declared here (rather than left to spring into existence inside the
% nested function) purely for readability.
selChanIdx = initChanIdx;    % indices into chLabelsFull
selICIdx   = [];             % component indices, into icaPrimary's decomposition
if ~isempty(icaPrimary)
    selICIdx = 1:min(10, src(icaPrimary).nComp);   % "top 10 ICs by default"
end
nRows      = 0;
spacing    = spacingChan;
offsets    = [];
yLimits    = [0 1];
ampUnit    = [char(956) 'V'];
zeroLines  = gobjects(0, 1);

%%% --- Figure ---
% Plain figure+axes (not uifigure/App Designer), 'painters' renderer, no
% antialiasing: on a machine with no GPU, OpenGL falls back to slow
% software emulation. (If you want to A/B test, 'opengl' is the only other
% option worth trying -- but painters is the safer default without a GPU.)
if nSrc == 2
    figName = sprintf('%s  vs  %s', src(1).name, src(2).name);
else
    figName = src(1).name;
end
fig = figure('Name', figName, 'NumberTitle', 'off', ...
    'Color', 'w', 'MenuBar', 'none', 'ToolBar', 'none', ...
    'Renderer', 'painters', 'GraphicsSmoothing', 'off', ...
    'Units', 'normalized', 'Position', [0.15 0.15 0.70 0.70], ...
    'WindowKeyPressFcn', @onKeyPress);

% panelX/panelW are shared by ax, axHypno, and axScroll so all three line
% up horizontally, starting near the actual left edge of the figure.
panelX = 0.05; panelW = 0.93;
% Vertical budget, bottom-to-top: control strip (0.012-0.077), scrollbar
% (0.089-0.119), hypnogram + its xlabel (0.150-0.205), then this panel --
% whose own tick labels and xlabel need the gap below 0.265 to itself.
ax = axes('Parent', fig, 'Units', 'normalized', 'Position', [panelX 0.265 panelW 0.700]);
ax.Box        = 'on';
ax.TickLength = [0 0];
ax.Toolbar    = [];
ax.SortMethod = 'childorder';   % skip depth-sorting -- irrelevant for 2D line data
disableDefaultInteractivity(ax);
ax.XTickMode  = 'auto';
% YTick/YTickLabel/YLim and the per-row zero-lines depend on the current
% selection, so they're (re)built by applyRowSelection() instead of being
% set once here.

% Vertical time-grid reference lines: one line object encoding multiple
% whole-second gridlines (NaN-separated, same trick as the multi-row signal
% lines below), full-height so they pass through every row's zero baseline.
% Transparent grey (RGBA colour) so they read as a subtle timing aid, not
% competing with the traces. Created once, XData/YData rebuilt every redraw
% since which seconds are in view changes with navigation. HandleVisibility
% off keeps it (and everything below) out of the legend -- see
% applyRowSelection()/redraw() for why that matters.
gridLine = line(ax, NaN, NaN, 'Color', [0.5 0.5 0.5 0.35], 'LineWidth', 0.5, ...
    'HitTest', 'off', 'PickableParts', 'none', 'HandleVisibility', 'off');

% Each trace is ONE line spanning all rows (NaN gaps lift the pen between
% rows) -- created once, updated via XData/YData, never deleted/recreated.
% This is what a cla()+plot() approach was paying for every single redraw.
% There is one line object per possible trace (raw + cleaned per source);
% unused ones are simply hidden.
traceLines = gobjects(maxSeries, 1);
for j = 1:maxSeries
    traceLines(j) = line(ax, NaN, NaN, 'Color', royalBlue, 'LineWidth', 0.5, ...
        'HitTest', 'off', 'PickableParts', 'none', 'Visible', 'off', ...
        'Tag', sprintf('eegTrace%d', j));
end

% Amplitude scale bar: a vertical reference line + label, updated (not
% recreated) every redraw since it depends on both ampScale and window.
% HandleVisibility off: legend AutoUpdate defaults to 'on', which would
% otherwise silently add "dataN" legend entries for these the moment
% they're (re)created after the legend already exists.
scaleBarLine = line(ax, [NaN NaN], [NaN NaN], 'Color', 'k', 'LineWidth', 1.2, ...
    'HitTest', 'off', 'PickableParts', 'none', 'HandleVisibility', 'off');
scaleBarText = text(ax, NaN, NaN, '', 'VerticalAlignment', 'middle', ...
    'FontSize', 8, 'Color', 'k', 'HitTest', 'off');

xlabel(ax, 'Time (s)');

%%% --- Controls: one strip of labelled groups along the bottom ---
% Related controls sit together in a titled uipanel instead of running
% along one long undifferentiated row. Two things keep it aligned:
%   * every panel's width is derived from the same weight table as its
%     contents, so one "unit" is the same number of pixels in EVERY panel
%     and buttons line up across group boundaries;
%   * each panel places its contents from a weight vector (layoutSlots)
%     rather than a hand-carried x-cursor, which is what let the old strip
%     drift out of true every time a control was added or relabelled.
% Order bottom-to-top: controls, scrollbar, hypnogram, EEG panel.
icaEnable = 'off';
if ~isempty(icaSrc), icaEnable = 'on'; end

ctrlY = 0.012; ctrlH = 0.065;
slotGap  = 0.10;                            % gap between slots, in slot-width units
navW     = [1.15 1.15 0.75 0.55 1.15 0.55];   % Prev  Next  "Epoch:" []  "Window (s):" []
stageW   = [1.30 0.80 0.80];                  % [stage]  First  Next
dispW    = [1.25 0.60 0.60 1.20];             % [time unit]  Amp-  Amp+  Channels...
traceW   = [1.00 1.30];                       % Diff  Clean only
icaW     = [1.15 1.30];                       % Show ICs  Subtract ICs
groupW   = cellfun(@(w) sum(w) + slotGap*(numel(w)-1), {navW, stageW, dispW, traceW, icaW});

panGap = 0.008; panX0 = 0.02;
groupW = groupW / sum(groupW) * ((0.98 - panX0) - panGap*(numel(groupW)-1));
groupX = panX0 + [0, cumsum(groupW(1:end-1) + panGap)];
yB = 0.10; hB = 0.74;    % buttons / popups inside a panel
yT = 0.04; hT = 0.62;    % text labels sit lower -- uicontrol text is top-aligned in its box

%--- Navigate ---
p = mkPanel(1, 'Navigate');
s = layoutSlots(navW, slotGap);
uicontrol(p, 'Style', 'pushbutton', 'String', char(9664) + " Prev", 'Units', 'normalized', ...
    'Position', [s(1,1) yB s(1,2) hB], 'Callback', @(o,e) shiftEpoch(-1));
uicontrol(p, 'Style', 'pushbutton', 'String', "Next " + char(9654), 'Units', 'normalized', ...
    'Position', [s(2,1) yB s(2,2) hB], 'Callback', @(o,e) shiftEpoch(+1));
uicontrol(p, 'Style', 'text', 'String', 'Epoch:', 'Units', 'normalized', ...
    'Position', [s(3,1) yT s(3,2) hT], 'HorizontalAlignment', 'right', 'BackgroundColor', 'w');
uicontrol(p, 'Style', 'edit', 'String', '1', 'Units', 'normalized', ...
    'Position', [s(4,1) yB s(4,2) hB], ...
    'Callback', @(o,e) jumpToEpoch(round(str2double(o.String))));
% NOTE: that edit box is input-only -- it is deliberately NOT synced back to
% the current epoch on every redraw (that's what the title is for). That
% sync was one of the uicontrol touches that turned out to be expensive.
uicontrol(p, 'Style', 'text', 'String', 'Window (s):', 'Units', 'normalized', ...
    'Position', [s(5,1) yT s(5,2) hT], 'HorizontalAlignment', 'right', 'BackgroundColor', 'w');
uicontrol(p, 'Style', 'edit', 'String', num2str(displaySec), 'Units', 'normalized', ...
    'Position', [s(6,1) yB s(6,2) hB], 'Callback', @(o,e) setDisplaySec(str2double(o.String)));

%--- Sleep stage ---
p = mkPanel(2, 'Sleep stage');
s = layoutSlots(stageW, slotGap);
popStage = uicontrol(p, 'Style', 'popupmenu', 'String', stageNames, 'Units', 'normalized', ...
    'Position', [s(1,1) yB s(1,2) hB]);
uicontrol(p, 'Style', 'pushbutton', 'String', 'First', 'Units', 'normalized', ...
    'Position', [s(2,1) yB s(2,2) hB], 'Callback', @(o,e) jumpToStage('first'));
uicontrol(p, 'Style', 'pushbutton', 'String', 'Next', 'Units', 'normalized', ...
    'Position', [s(3,1) yB s(3,2) hB], 'Callback', @(o,e) jumpToStage('next'));

%--- Display ---
p = mkPanel(3, 'Display');
s = layoutSlots(dispW, slotGap);
popScale = uicontrol(p, 'Style', 'popupmenu', 'String', {'seconds', 'minutes', 'hours'}, ...
    'Units', 'normalized', 'Position', [s(1,1) yB s(1,2) hB], 'Callback', @(o,e) redraw());
uicontrol(p, 'Style', 'pushbutton', 'String', 'Amp -', 'Units', 'normalized', ...
    'Position', [s(2,1) yB s(2,2) hB], 'Callback', @(o,e) scaleDown());
uicontrol(p, 'Style', 'pushbutton', 'String', 'Amp +', 'Units', 'normalized', ...
    'Position', [s(3,1) yB s(3,2) hB], 'Callback', @(o,e) scaleUp());
btnRows = uicontrol(p, 'Style', 'pushbutton', 'String', 'Channels...', 'Units', 'normalized', ...
    'Position', [s(4,1) yB s(4,2) hB], 'Callback', @(o,e) openRowSelector());

%--- Traces ---
p = mkPanel(4, 'Traces');
s = layoutSlots(traceW, slotGap);
btnDiff = uicontrol(p, 'Style', 'togglebutton', 'String', 'Diff (d)', 'Units', 'normalized', ...
    'Position', [s(1,1) yB s(1,2) hB], 'Callback', @(o,e) toggleDiff());
btnClean = uicontrol(p, 'Style', 'togglebutton', 'String', 'Clean only (c)', 'Units', 'normalized', ...
    'Position', [s(2,1) yB s(2,2) hB], 'Callback', @(o,e) toggleCleanOnly());

%--- ICA ---
p = mkPanel(5, 'ICA');
s = layoutSlots(icaW, slotGap);
% The label names the ACTION, not the current state, so it flips to
% "Show chans" once component activations are on screen (see toggleICView).
btnICView = uicontrol(p, 'Style', 'togglebutton', 'String', 'Show ICs (i)', 'Units', 'normalized', ...
    'Position', [s(1,1) yB s(1,2) hB], 'Enable', icaEnable, 'Callback', @(o,e) toggleICView());
btnSubtract = uicontrol(p, 'Style', 'togglebutton', 'String', 'Subtract ICs', 'Units', 'normalized', ...
    'Position', [s(2,1) yB s(2,2) hB], 'Enable', icaEnable, 'Callback', @(o,e) toggleSubtract());

%%% --- Hypnogram strip (static, drawn once -- never touched in redraw) ---
% Shares panelX/panelW with the EEG panel and the scrollbar, so all three
% line up horizontally.
hypnoLabels = {'N3', 'N2', 'N1', 'R', 'W'};   % short form of stageNames, same order as stageValues
hypnoY = 0.150; hypnoH = 0.055;
axHypno = axes('Parent', fig, 'Units', 'normalized', ...
    'Position', [panelX, hypnoY, panelW, hypnoH]);

if ~isempty(scoringDigits)
    % XLim ends exactly at the last scored epoch -- not at totalDur, which
    % is the EEG recording's own length and isn't guaranteed to match the
    % scored range (that mismatch was leaving empty space at the end).
    hypnoEnd   = numel(scoringDigits) * scoreEpochSec;
    epochTimes = (0:numel(scoringDigits)-1) * scoreEpochSec;
    stairs(axHypno, [epochTimes, hypnoEnd] / 3600, double([scoringDigits(:); scoringDigits(end)]), ...
        'Color', [0.20 0.20 0.55], 'LineWidth', 1);
    hypnoXLim = [0, hypnoEnd] / 3600;
else
    hypnoXLim = [0, max(totalDur, displaySec)] / 3600;
end

% Axes cosmetics go AFTER the stairs call: plotting into an axes resets its
% tick labels to auto, which was quietly turning the stage names back into
% the raw scoring digits (-3..1) on the y-axis.
axHypno.YLim       = [min(stageValues)-0.5, max(stageValues)+0.5];
axHypno.YTick      = stageValues;
axHypno.YTickLabel = hypnoLabels;
axHypno.XLim       = hypnoXLim;
axHypno.XTickMode  = 'auto';
axHypno.Box        = 'on';
axHypno.FontSize   = 7;
axHypno.TickLength = [0 0];
axHypno.Toolbar    = [];
disableDefaultInteractivity(axHypno);
xlabel(axHypno, 'Time (h)', 'FontSize', 7);

if isempty(scoringDigits)
    text(axHypno, mean(axHypno.XLim), mean(axHypno.YLim), 'No sleep scoring loaded', ...
        'HorizontalAlignment', 'center', 'FontSize', 8, 'Color', [0.5 0.5 0.5]);
end

%%% --- Custom scrollbar ---
% Drawn as patches inside a thin axes instead of a native uicontrol slider.
% On a GPU-less machine, moving a native Win32/Java slider turned out to
% cost 700-1300ms per redraw (confirmed by timing it in isolation from
% everything else), dwarfing the ~40ms the actual line drawing takes. A
% patch is rendered by the same painters pipeline as the traces and costs
% next to nothing by comparison.
sliderY = 0.089; sliderH = 0.030;   % clears the control strip below
axScroll = axes('Parent', fig, 'Units', 'normalized', ...
    'Position', [panelX, sliderY, panelW, sliderH]);
axScroll.XLim       = [0, max(totalDur, displaySec)];
axScroll.YLim       = [0 1];
axScroll.XTick      = [];
axScroll.YTick      = [];
axScroll.Box        = 'on';
axScroll.Color      = [0.90 0.90 0.90];
axScroll.Toolbar    = [];
disableDefaultInteractivity(axScroll);

patch(axScroll, [0 totalDur totalDur 0], [0 0 1 1], [0.90 0.90 0.90], ...
    'EdgeColor', 'none', 'ButtonDownFcn', @scrollDown);
scrollThumb = patch(axScroll, [0 displaySec displaySec 0], [0 0 1 1], [0.55 0.65 0.88], ...
    'EdgeColor', [0.30 0.35 0.55], 'LineWidth', 1, 'ButtonDownFcn', @scrollDown);

applyRowSelection();   % builds ticks/ylim/zero-lines and redraws


%% ---- nested callbacks (share this function's workspace directly) ----

    function p = mkPanel(i, ttl)
        % One titled group in the bottom control strip. Width comes from
        % groupW, which is derived from the same weight table the panel's
        % contents use -- that shared scale is what keeps buttons aligned
        % across panel boundaries.
        p = uipanel(fig, 'Title', ttl, 'Units', 'normalized', ...
            'Position', [groupX(i) ctrlY groupW(i) ctrlH], ...
            'BackgroundColor', 'w', 'FontSize', 8, 'TitlePosition', 'centertop');
    end

    function onKeyPress(~, evt)
        tCb = tic;
        switch evt.Key
            case 'leftarrow',  shiftEpoch(-1);
            case 'rightarrow', shiftEpoch(+1);
            case 'd',          toggleDiff();
            case 'c',          toggleCleanOnly();
            case 'i',          toggleICView();
            otherwise, return
        end
        fprintf('[onKeyPress:%s] TOTAL %6.1f ms\n\n', evt.Key, toc(tCb)*1000);
    end

    function scrollDown(~, ~)
        % A single discrete jump per click -- deliberately NOT wired to
        % WindowButtonMotionFcn, so the epoch only changes on an actual
        % click, not on every mouse move while the button happens to be
        % held down or the cursor passes over the scrollbar.
        tCb = tic;
        cp     = fig.CurrentPoint;
        axPos  = axScroll.Position;
        frac   = max(0, min(1, (cp(1) - axPos(1)) / axPos(3)));
        clickT = frac * axScroll.XLim(2);
        windowStart = clickT - displaySec/2;
        redraw();
        fprintf('[scrollDown] TOTAL %6.1f ms\n\n', toc(tCb)*1000);
    end

    function shiftEpoch(delta)
        curEpoch = floor(windowStart / scoreEpochSec) + 1;
        jumpToEpoch(curEpoch + delta);
    end

    function jumpToEpoch(epochIdx)
        tCb = tic;
        if ~isfinite(epochIdx), return; end
        epochIdx = max(1, min(nEpochs, epochIdx));
        windowStart = (epochIdx - 1) * scoreEpochSec;
        redraw();
        fprintf('[jumpToEpoch] TOTAL %6.1f ms\n\n', toc(tCb)*1000);
    end

    function toggleDiff()
        tCb = tic;
        diffMode = ~diffMode;
        if diffMode, cleanOnly = false; end   % diff and clean-only are mutually exclusive display modes
        redraw();
        fprintf('[toggleDiff] TOTAL %6.1f ms\n\n', toc(tCb)*1000);
    end

    function toggleCleanOnly()
        tCb = tic;
        cleanOnly = ~cleanOnly;
        if cleanOnly, diffMode = false; end   % see toggleDiff()
        redraw();
        fprintf('[toggleCleanOnly] TOTAL %6.1f ms\n\n', toc(tCb)*1000);
    end

    function toggleICView()
        % Channel traces <-> component activations. Inert without an ICA
        % decomposition (the button is greyed out too, this guards the
        % keyboard shortcut).
        if isempty(icaPrimary), return; end
        tCb = tic;
        if strcmp(viewMode, 'chan')
            viewMode = 'ic';
            set(btnRows,   'String', 'ICs...');
            set(btnICView, 'String', 'Show chans (i)');
        else
            viewMode = 'chan';
            set(btnRows,   'String', 'Channels...');
            set(btnICView, 'String', 'Show ICs (i)');
        end
        applyRowSelection();
        fprintf('[toggleICView:%s] TOTAL %6.1f ms\n\n', viewMode, toc(tCb)*1000);
    end

    function toggleSubtract()
        % Switching subtraction ON opens the component picker first, so the
        % set being projected out is always something the user just
        % confirmed. Switching it OFF just hides the extra trace -- the
        % selection is remembered, so toggling back on and confirming
        % re-applies it. (To EDIT the set while it is on: toggle off, then
        % on again.)
        if isempty(icaPrimary), set(btnSubtract, 'Value', 0); return; end
        tCb = tic;
        if subtractOn
            subtractOn = false;
        elseif openSubtractSelector()
            subtractOn = any(arrayfun(@(k) ~isempty(src(k).subIdx), icaSrc));
        end
        set(btnSubtract, 'Value', subtractOn);
        redraw();
        fprintf('[toggleSubtract] TOTAL %6.1f ms\n\n', toc(tCb)*1000);
    end

    function confirmed = openSubtractSelector()
        % One picker per ICA-bearing source (they are separate
        % decompositions, so they need separate component sets).
        % NOTE on loop variables in this file's nested functions: they share
        % the parent workspace, so a bare `for k` here would clobber the
        % setup loops' k. Hence kSub/kSer/jLine below.
        confirmed = false;
        for kSub = icaSrc
            [sel, ok] = viewer.pickerDialog(src(kSub).icText, src(kSub).subIdx, ...
                sprintf('Subtract ICs - %s', src(kSub).name), ...
                'Components to project out (highlighted = currently subtracted):', ...
                'Subtract', src(kSub).classNames, src(kSub).icClass);
            if ~ok, continue; end
            src(kSub).subIdx = sel;
            % Fold the whole "unmix, zero the bad components, remix" round
            % trip into ONE channels-by-channels operator, so redraw() only
            % ever does a small matmul on the visible window instead of
            % materialising full-length activations or a full cleaned copy
            % of the recording.
            src(kSub).subP = src(kSub).winv(:, sel) * src(kSub).unmix(sel, :);
            confirmed = true;
        end
    end

    function openRowSelector()
        if strcmp(viewMode, 'chan')
            [sel, ok] = viewer.pickerDialog(chLabelsFull, selChanIdx, 'Select channels', ...
                'Channels to plot (highlighted = currently shown):', ...
                'Plot selected', {}, []);
            if ~ok || isempty(sel), return; end   % cancelled, or selection cleared -- leave plot as-is
            selChanIdx = sel;
        else
            p = icaPrimary;
            [sel, ok] = viewer.pickerDialog(src(p).icText, selICIdx, 'Select components', ...
                'Components to plot (highlighted = currently shown):', ...
                'Plot selected', src(p).classNames, src(p).icClass);
            if ~ok || isempty(sel), return; end
            selICIdx = sel;
        end
        applyRowSelection();
    end

    function applyRowSelection()
        % Rebuilds everything that depends on WHICH rows are stacked in the
        % panel: y-ticks, y-limits, per-row zero-lines. Called on startup,
        % on a view-mode switch, and on a channel/IC selection change.
        tCb = tic;
        if strcmp(viewMode, 'chan')
            rowLabels = chLabelsFull(selChanIdx);
            spacing   = spacingChan;
            ampUnit   = [char(956) 'V'];
        else
            rowLabels = src(icaPrimary).icTick(selICIdx);
            spacing   = spacingIC;
            ampUnit   = 'a.u.';
        end
        nRows = numel(rowLabels);
        if nRows == 0
            fprintf('eegCompareViewer: row selection cannot be empty; ignoring.\n');
            return
        end

        offsets = (nRows - 1 : -1 : 0)' * spacing;   % first row = top
        yLimits = [-spacing*0.75, offsets(1) + spacing*0.75];

        ax.YTick      = flip(offsets);
        ax.YTickLabel = flip(rowLabels);
        ax.YLim       = yLimits;

        % Row count changed -- the per-row zero-lines have to be rebuilt
        % from scratch (the trace lines are unaffected: they're single line
        % objects whose XData/YData redraw() rewrites wholesale).
        delete(zeroLines(isgraphics(zeroLines)));
        zeroLines = gobjects(nRows, 1);
        for iRow = 1:nRows
            % HandleVisibility off: without it, legend AutoUpdate ('on' by
            % default) treats each freshly-created line as a new plot to
            % auto-label -- that's what was showing up as "data1"/"data2"
            % entries in the legend every time the selection changed.
            zeroLines(iRow) = line(ax, [-1e6 1e6], [offsets(iRow) offsets(iRow)], ...
                'Color', zeroLineGray, 'LineStyle', ':', 'LineWidth', 0.5, ...
                'HitTest', 'off', 'PickableParts', 'none', 'HandleVisibility', 'off');
        end
        uistack(traceLines, 'top');   % keep the traces above the freshly-added zero-lines

        redraw();
        fprintf('[applyRowSelection] TOTAL %6.1f ms\n\n', toc(tCb)*1000);
    end

    function scaleUp()
        ampScale = min(ampScale * 1.25, 40);
        redraw();
    end

    function scaleDown()
        ampScale = max(ampScale / 1.25, 0.025);
        redraw();
    end

    function setDisplaySec(val)
        if ~isfinite(val) || val <= 0, return; end
        displaySec = min(val, max(totalDur, 1));
        axScroll.XLim = [0, max(totalDur, displaySec)];
        if isempty(scoringDigits)
            % Only the no-scoring placeholder range depends on displaySec;
            % with real scoring loaded, the hypnogram's span is fixed to
            % the scored range regardless of the current window size.
            axHypno.XLim = [0, max(totalDur, displaySec)] / 3600;
        end
        set(scrollThumb, 'XData', [0, displaySec, displaySec, 0]);
        redraw();
    end

    function jumpToStage(mode)
        tCb = tic;
        if isempty(scoringDigits)
            fprintf('eegCompareViewer: no scoring file loaded.\n');
            return
        end
        stageIdx = get(popStage, 'Value');
        val      = stageValues(stageIdx);
        curEpoch = floor(windowStart / scoreEpochSec) + 1;

        if strcmp(mode, 'first')
            idx = find(scoringDigits == val, 1, 'first');
        else
            idx = find(scoringDigits == val & (1:numel(scoringDigits)) > curEpoch, 1, 'first');
        end

        if isempty(idx)
            fprintf('eegCompareViewer: no further "%s" epoch found.\n', stageNames{stageIdx});
            return
        end
        jumpToEpoch(idx);
        fprintf('[jumpToStage:%s] TOTAL %6.1f ms\n\n', mode, toc(tCb)*1000);
    end

    function redraw()
        tTotal = tic;

        tStep = tic;
        windowStart = max(0, min(windowStart, max(totalDur - displaySec, 0)));
        idxStart = max(1, round(windowStart * srate) + 1);
        idxEnd   = min(nSamp, idxStart + round(displaySec * srate) - 1);
        t        = (idxStart:idxEnd) / srate;
        % Only every plotDecimation-th sample is actually drawn -- this is
        % display-only downsampling (xlim/title/scale bar still use the
        % full-resolution t/tScaled), it doesn't touch src(k).data.
        plotIdx  = idxStart:plotDecimation:idxEnd;
        fprintf('  [slice]          %6.1f ms\n', toc(tStep)*1000);

        % Time-unit conversion: divide the data itself into the selected
        % unit and let MATLAB auto-generate tick labels natively, instead
        % of hand-formatting every tick label every redraw.
        scaleOpts = {'seconds', 'minutes', 'hours'};
        switch scaleOpts{get(popScale, 'Value')}
            case 'seconds', divisor = 1;    unitStr = 's';
            case 'minutes', divisor = 60;   unitStr = 'min';
            case 'hours',   divisor = 3600; unitStr = 'h';
        end
        tScaled     = t / divisor;
        tScaledPlot = tScaled(1:plotDecimation:end);

        %%% --- Assemble the traces for this window ---
        % Order matters: raw then cleaned, source 1 then source 2. That
        % makes series(1)-series(2) the meaningful difference in every
        % configuration (two inputs -> input1 vs input2; one input with
        % subtraction -> what the components carried away), and makes the
        % LAST series the cleanest one, which is what "Clean only" shows.
        tStep = tic;
        series = struct('y', {}, 'color', {}, 'name', {});
        if strcmp(viewMode, 'chan')
            for kSer = 1:nSrc
                rows = src(kSer).rowOf(selChanIdx);
                raw  = src(kSer).data(rows, plotIdx);
                series(end+1) = struct('y', raw, 'color', colRaw{kSer}, ...
                    'name', src(kSer).name); %#ok<AGROW>
                if subtractOn && src(kSer).hasICA && ~isempty(src(kSer).subIdx)
                    % Project the selected components out, for the visible
                    % rows and window only -- see openSubtractSelector()
                    % for why subP exists.
                    pos = src(kSer).posInICA(selChanIdx);
                    hit = pos > 0;
                    cleaned = raw;
                    if any(hit)
                        cleaned(hit, :) = raw(hit, :) - ...
                            src(kSer).subP(pos(hit), :) * src(kSer).data(src(kSer).icachansind, plotIdx);
                    end
                    series(end+1) = struct('y', cleaned, 'color', colClean{kSer}, ...
                        'name', sprintf('%s  (-%d IC)', src(kSer).name, numel(src(kSer).subIdx))); %#ok<AGROW>
                end
            end
        else
            for kSer = icaSrc
                if max(selICIdx) > src(kSer).nComp
                    continue   % this decomposition has fewer components than the selection
                end
                act = src(kSer).unmix(selICIdx, :) * src(kSer).data(src(kSer).icachansind, plotIdx);
                series(end+1) = struct('y', act, 'color', colRaw{kSer}, ...
                    'name', sprintf('%s  IC activations', src(kSer).name)); %#ok<AGROW>
            end
        end

        if diffMode && numel(series) >= 2
            series = struct('y', series(1).y - series(2).y, 'color', diffColor, ...
                'name', sprintf('Difference (%s - %s)', series(1).name, series(2).name));
        elseif cleanOnly && numel(series) >= 2
            series = series(end);
        end
        fprintf('  [series]         %6.1f ms\n', toc(tStep)*1000);

        tStep = tic;
        Xall = repmat([tScaledPlot, NaN], 1, nRows);
        nS   = numel(series);
        for jLine = 1:maxSeries
            if jLine <= nS
                Yall = reshape([series(jLine).y * ampScale + offsets, nan(nRows, 1)].', 1, []);
                set(traceLines(jLine), 'XData', Xall, 'YData', Yall, ...
                    'Color', series(jLine).color, 'Visible', 'on');
            else
                set(traceLines(jLine), 'Visible', 'off');
            end
        end
        fprintf('  [lines]          %6.1f ms\n', toc(tStep)*1000);

        tStep = tic;
        % Vertical whole-second reference gridlines across the current
        % window, full axes height (so they run through every row's zero
        % baseline). Step size auto-coarsens for big windows so this never
        % tries to paint hundreds of lines (1 s step for the default 30 s
        % epoch view, matching "a line at every displayed second").
        gStep = gridStepCandidates(find(displaySec ./ gridStepCandidates <= maxGridLines, 1, 'first'));
        if isempty(gStep), gStep = gridStepCandidates(end); end
        gridT = ((ceil(t(1) / gStep) * gStep) : gStep : t(end)) / divisor;
        nGrid = numel(gridT);
        set(gridLine, 'XData', reshape([gridT; gridT; nan(1, nGrid)], 1, []), ...
                      'YData', repmat([yLimits(1), yLimits(2), NaN], 1, nGrid));
        fprintf('  [time grid]      %6.1f ms\n', toc(tStep)*1000);

        tStep = tic;
        % Scale bar: a "nice" 1/2/5-times-a-power-of-ten value near 30% of
        % the row spacing. Derived rather than picked from a fixed list,
        % because IC activations live nowhere near the microvolt range.
        barVal    = niceStep(0.3 * spacing / ampScale);
        barHeight = barVal * ampScale;
        barX      = tScaled(1) + 0.015 * (tScaled(end) - tScaled(1));
        barY0     = yLimits(2) - barHeight - spacing * 0.05;
        set(scaleBarLine, 'XData', [barX barX], 'YData', [barY0, barY0 + barHeight]);
        set(scaleBarText, 'Position', [barX + 0.01*(tScaled(end)-tScaled(1)), barY0 + barHeight/2, 0], ...
            'String', sprintf('%g %s', barVal, ampUnit));
        fprintf('  [scale bar]      %6.1f ms\n', toc(tStep)*1000);

        tStep = tic;
        xlim(ax, [tScaled(1), tScaled(end)]);
        xlabel(ax, sprintf('Time (%s)', unitStr));
        fprintf('  [xlim]           %6.1f ms\n', toc(tStep)*1000);

        %%% --- Title: epoch, time range, sleep stage, mode ---
        tStep = tic;
        curEpoch = floor(windowStart / scoreEpochSec) + 1;
        stageStr = '';
        if ~isempty(scoringDigits) && curEpoch <= numel(scoringDigits)
            hit = find(stageValues == scoringDigits(curEpoch), 1);
            if ~isempty(hit)
                stageStr = sprintf('   |   Stage: %s', stageNames{hit});
            end
        end
        modeStr = '';
        if strcmp(viewMode, 'ic'), modeStr = [modeStr '   |   IC ACTIVATIONS']; end
        if diffMode
            modeStr = [modeStr '   |   DIFFERENCE'];
        elseif cleanOnly && nS >= 1
            modeStr = [modeStr sprintf('   |   ONLY: %s', series(end).name)];
        end
        title(ax, sprintf('Epoch %d / %d   (%.1f-%.1f s)%s%s', ...
            curEpoch, nEpochs, t(1), t(end), stageStr, modeStr), 'Interpreter', 'none');
        fprintf('  [title]          %6.1f ms\n', toc(tStep)*1000);

        % Legend + mode buttons: only touched when the display state
        % actually changes, not on every pure-navigation redraw (native
        % uicontrol/legend churn is the expensive part on a GPU-less box).
        % AutoUpdate off so MATLAB never auto-adds a "dataN" entry for any
        % object created after this call (e.g. the per-row zero-lines
        % rebuilt on a selection change) -- HandleVisibility off on those
        % already prevents it too, but belt-and-suspenders costs nothing.
        tStep = tic;
        stateKey = sprintf('%s|%d%d%d|%s', viewMode, diffMode, cleanOnly, subtractOn, ...
            strjoin({series.name}, '|'));
        if ~strcmp(stateKey, lastStateKey)
            if nS > 0
                legend(ax, traceLines(1:nS), {series.name}, 'Location', 'northeast', ...
                    'Interpreter', 'none', 'AutoUpdate', 'off');
            else
                legend(ax, 'off');
            end
            set(btnDiff,     'Value', diffMode);
            set(btnClean,    'Value', cleanOnly);
            set(btnSubtract, 'Value', subtractOn);
            set(btnICView,   'Value', strcmp(viewMode, 'ic'));
            lastStateKey = stateKey;
        end
        fprintf('  [legend/buttons] %6.1f ms\n', toc(tStep)*1000);

        tStep = tic;
        set(scrollThumb, 'XData', [windowStart, windowStart+displaySec, windowStart+displaySec, windowStart]);
        fprintf('  [scrollbar]      %6.1f ms\n', toc(tStep)*1000);

        tStep = tic;
        drawnow
        fprintf('  [drawnow]        %6.1f ms\n', toc(tStep)*1000);

        fprintf('  TOTAL redraw     %6.1f ms\n\n', toc(tTotal)*1000);
    end

end


%% ---- local helpers (no access to the viewer's workspace) ----

function s = attachICA(s, EEG)
% Pull the ICA decomposition (and, if present, the ICLabel classification)
% out of an EEGLAB structure into the viewer's own per-source record.
% Everything stays in single precision and NOTHING full-length is
% materialised here: activations and cleaned data are both derived per
% display window inside redraw().
s.hasICA     = false;
s.nComp      = 0;
s.icTick     = {};
s.icText     = {};
s.icClass    = [];
s.classNames = {};
s.subIdx     = [];

if ~isfield(EEG, 'icaweights') || isempty(EEG.icaweights) || ...
   ~isfield(EEG, 'icasphere')  || isempty(EEG.icasphere)
    return
end

if isfield(EEG, 'icachansind') && ~isempty(EEG.icachansind)
    icachansind = EEG.icachansind(:).';
else
    icachansind = 1:size(EEG.icasphere, 2);
end
% chans1020/pop_select upstream can reorder or drop channels without
% fixing icachansind. Refusing the decomposition is the honest response --
% silently indexing the wrong rows would produce plausible-looking traces
% of the wrong thing.
if max(icachansind) > size(EEG.data, 1) || numel(icachansind) ~= size(EEG.icasphere, 2)
    warning('eegCompareViewer:staleICA', ...
        'icachansind does not match the channel set (%d channels, %d sphere columns, max index %d); ignoring the ICA decomposition for "%s".', ...
        size(EEG.data, 1), size(EEG.icasphere, 2), max(icachansind), s.name);
    return
end

s.unmix       = single(EEG.icaweights * EEG.icasphere);   % nComp x nICAchan
s.icachansind = icachansind;
s.nComp       = size(s.unmix, 1);
if isfield(EEG, 'icawinv') && ~isempty(EEG.icawinv)
    s.winv = single(EEG.icawinv);                         % nICAchan x nComp
else
    s.winv = single(pinv(double(s.unmix)));
end

% Position of every common channel within the decomposition (0 = that
% channel was not part of the ICA and so cannot be cleaned).
s.posInICA = zeros(1, numel(s.rowOf));
[tf, loc]  = ismember(s.rowOf, icachansind);
s.posInICA(tf) = loc(tf);

% Component labels: "IC 12   Muscle  87%" in the pickers, "IC 12 (Muscle)"
% on the y-axis where space is tight.
s.icTick = arrayfun(@(i) sprintf('IC %d', i), 1:s.nComp, 'UniformOutput', false);
s.icText = s.icTick;
s.icClass = zeros(1, s.nComp);
if isfield(EEG, 'etc') && isfield(EEG.etc, 'ic_classification') && ...
        isfield(EEG.etc.ic_classification, 'ICLabel')
    L = EEG.etc.ic_classification.ICLabel;
    s.classNames = L.classes(:).';
    [p, cls] = max(L.classifications, [], 2);
    n = min(s.nComp, numel(cls));
    s.icClass(1:n) = cls(1:n);
    for i = 1:n
        s.icTick{i} = sprintf('IC %d (%s)', i, s.classNames{cls(i)});
        s.icText{i} = sprintf('IC %-4d %-11s %3.0f%%', i, s.classNames{cls(i)}, p(i)*100);
    end
end

% Default subtraction set: every non-brain component when ICLabel says so,
% otherwise whatever a previous selectcomps/pop_selectcomps run flagged.
if ~isempty(s.classNames)
    brainCls = find(strcmpi(s.classNames, 'Brain'), 1);
    s.subIdx = find(s.icClass > 0 & s.icClass ~= brainCls);
elseif isfield(EEG, 'reject') && isfield(EEG.reject, 'gcompreject') && ~isempty(EEG.reject.gcompreject)
    s.subIdx = find(EEG.reject.gcompreject(:).');
end
s.subP   = s.winv(:, s.subIdx) * s.unmix(s.subIdx, :);
s.hasICA = true;
end


function s = layoutSlots(weights, gap)
% Turn a vector of relative widths into normalized [x, width] rows inside a
% uipanel, with `gap` (in the same units as the weights) between each pair
% and a small inset at both ends. Every control strip in this file is
% positioned through here, so the slot geometry is defined once instead of
% being re-derived by hand at each call site.
inset = 0.012;
total = sum(weights) + gap * (numel(weights) - 1);
x     = [0, cumsum(weights(1:end-1) + gap)];
s     = [inset + x(:)/total * (1 - 2*inset), weights(:)/total * (1 - 2*inset)];
end


function v = robustSpacing(probes)
% Vertical gap between stacked rows: six robust standard deviations, so a
% single high-variance row can't squash everything else flat. Takes a cell
% of probe matrices and pools their per-row spreads -- the sources can have
% different row counts (two ICA decompositions need not have the same
% number of components), so they can't simply be concatenated.
sd = [];
for i = 1:numel(probes)
    sd = [sd; std(probes{i}, 0, 2)]; %#ok<AGROW>
end
v = 6 * median(sd, 'omitnan');
if ~isfinite(v) || v <= 0, v = 50; end
end


function v = niceStep(x)
% Nearest 1/2/5 x 10^n value to x -- used for the amplitude scale bar,
% which has to read sensibly in microvolts AND in the arbitrary units of
% component activations.
if ~isfinite(x) || x <= 0, v = 1; return; end
e = floor(log10(x));
m = x / 10^e;
if     m < 1.5, m = 1;
elseif m < 3.5, m = 2;
elseif m < 7.5, m = 5;
else,           m = 10;
end
v = m * 10^e;
end
