%% EEGCompareViewer
% Overlay two EEG recordings channel-by-channel for direct visual
% comparison, with epoch/stage navigation driven by a sleep scoring file.
% Requires EEGLAB (pop_interp, pop_select, ...) already on the path.
clc;
run('dependancies.m')

%% ===================== USER SETTINGS =====================
% file1       = '\\vs03.herseninstituut.knaw.nl\VS03-SandC-1\data\nin\data-drop\rawdata\derivatives\prep-ged\sub-drop0001\ses-t1\sub-drop0001_ses-t1_task-sleep_run-01_desc-zc2gedWakeBBAuto_eeg.set';   % first EEG file
file1       = '\\vs03.herseninstituut.knaw.nl\VS03-SandC-1\data\nin\data-drop\rawdata\derivatives\prep-ged\sub-drop0001\ses-t1\sub-drop0001_ses-t1_task-sleep_run-01_desc-zc_eeg.vhdr';   % first EEG file
file2       = '\\vs03.herseninstituut.knaw.nl\VS03-SandC-1\data\nin\data-drop\rawdata\derivatives\prep-ged\sub-drop0001\ses-t1\sub-drop0001_ses-t1_task-sleep_run-01_desc-zc2gedWakeBBAutoFSAutoPlus_eeg.set';   % second EEG file
% file2       = '\\vs03.herseninstituut.knaw.nl\VS03-SandC-1\data\nin\data-drop\rawdata\derivatives\prep-ged\sub-drop0001\ses-t1\sub-drop0001_ses-t1_task-sleep_run-01_desc-zc2gedWakeBBAutoPlusFSAutoPlus_eeg.set';   % second EEG file
scoringFile = '\\vs03.herseninstituut.knaw.nl\VS03-SandC-1\data\nin\data-drop\derivatives\scores\final\sub-drop0001_ses-t1_task-sleep_run-01_eeg.csv';                         % sleep scoring file, '' if none

net = 'EGI256'; % net type, passed to chans1020
targetrate = 100;

%% ===================== IMPORT + PREP =====================
EEG1 = fast_eeg_import(file1);
EEG2 = fast_eeg_import(file2);

EEG1 = prepEEG1020(EEG1, net, targetrate);
EEG2 = prepEEG1020(EEG2, net, targetrate);

scoringDigits = [];
if ~isempty(scoringFile)
    scoringDigits = scoreloader(scoringFile);
end

EEG1 = chans1020(EEG1, true, 'add_eog', true, 'net', net);
EEG2 = chans1020(EEG2, true, 'add_eog', true, 'net', net);

eegCompareViewer(EEG1, EEG2, scoringDigits, file1, file2);


%% ===================== LOCAL FUNCTIONS =====================

function EEG = prepEEG1020(EEG, net, targetrate)
    %%% Interpolate missing channels, then keep only 10-20 electrodes
    if targetrate < EEG.srate; EEG = pop_resample(EEG, targetrate); end
    try; EEG = pop_interp(EEG, EEG.urchanlocs, 'spherical'); end
%     EEG = chans1020(EEG, true, 'add_eog', true, 'net', net);
end


function eegCompareViewer(EEG1, EEG2, scoringDigits, file1, file2)
% All callbacks below are NESTED functions: they share this function's
% workspace directly (windowStart, diffMode, EEG1data, ...) instead of
% round-tripping state through get/setappdata on every redraw.

%%% --- Align channels between the two files (label match, stable order) ---
labels1 = {EEG1.chanlocs.labels};
labels2 = {EEG2.chanlocs.labels};
[~, ia, ib] = intersect(labels1, labels2, 'stable');
if numel(ia) < numel(labels1) || numel(ib) < numel(labels2)
    warning('eegCompareViewer:chanMismatch', ...
        'Channel labels differ between files; using %d common channels.', numel(ia));
end
EEG1.data     = EEG1.data(ia, :);
EEG1.chanlocs = EEG1.chanlocs(ia);
EEG2.data     = EEG2.data(ib, :);

if EEG1.srate ~= EEG2.srate
    warning('eegCompareViewer:srateMismatch', ...
        'Sampling rates differ (%.3f vs %.3f Hz); using file 1''s rate.', EEG1.srate, EEG2.srate);
end

chLabels = {EEG1.chanlocs.labels};
nCh      = numel(chLabels);
srate    = EEG1.srate;
nSamp    = min(size(EEG1.data, 2), size(EEG2.data, 2));
EEG1data = single(EEG1.data(:, 1:nSamp));
EEG2data = single(EEG2.data(:, 1:nSamp));
totalDur      = nSamp / srate;         % seconds
scoreEpochSec = 30;                    % fixed: sleep-scoring epoch length, used for epoch/stage alignment
nEpochs       = max(1, floor(totalDur / scoreEpochSec));

%%% --- Colours ---
royalBlue    = [ 65 105 225] / 255;
redColor     = [220  20  60] / 255;    % crimson -- a rich, readable red
diffColor    = [230 159   0] / 255;    % amber, distinct from both signal colours
zeroLineGray = [0.75 0.75 0.75];

%%% --- Amplitude scaling for stacked channel offsets (robust, subsampled) ---
probe   = [EEG1data(:, 1:50:end), EEG2data(:, 1:50:end)];
spacing = 6 * median(std(probe, 0, 2), 'omitnan');
if ~isfinite(spacing) || spacing <= 0, spacing = 50; end
offsets = (nCh - 1 : -1 : 0)' * spacing;   % channel 1 = top row
tickPos    = flip(offsets);
tickLabels = flip(chLabels);
yLimits    = [-spacing*0.75, offsets(1) + spacing*0.75];

%%% --- Sleep-stage lookup (qol/scoreloader.m convention) ---
stageNames  = {'N3', 'N2', 'N1', 'REM', 'Wake'};
stageValues = [-3, -2, -1, 0, 1];

%%% --- Navigation / display state (shared with all nested callbacks) ---
windowStart     = 0;
displaySec      = scoreEpochSec;   % how much time is shown at once; user-adjustable, defaults to 30 s
diffMode        = false;
lastDiffMode    = NaN;      % forces the legend/diff-button to be set once on first redraw
ampScale        = 1;        % vertical gain applied to the traces before offsetting
scaleCandidates = [2000 1000 500 200 100 50 20 10 5 2 1 0.5 0.2 0.1];   % "nice" uV steps for the scale bar

%%% --- Figure ---
% Plain figure+axes (not uifigure/App Designer), 'painters' renderer, no
% antialiasing: on a machine with no GPU, OpenGL falls back to slow
% software emulation. (If you want to A/B test, 'opengl' is the only other
% option worth trying -- but painters is the safer default without a GPU.)
[~, n1] = fileparts(file1); [~, n2] = fileparts(file2);
fig = figure('Name', sprintf('%s  vs  %s', n1, n2), 'NumberTitle', 'off', ...
    'Color', 'w', 'MenuBar', 'none', 'ToolBar', 'none', ...
    'Renderer', 'painters', 'GraphicsSmoothing', 'off', ...
    'Units', 'normalized', 'Position', [0.15 0.15 0.70 0.70], ...
    'WindowKeyPressFcn', @onKeyPress);

% panelX/panelW are shared by ax, axHypno, and axScroll so all three line
% up horizontally, starting near the actual left edge of the figure.
panelX = 0.05; panelW = 0.93;
ax = axes('Parent', fig, 'Units', 'normalized', 'Position', [panelX 0.20 panelW 0.77]);
ax.Box        = 'on';
ax.TickLength = [0 0];
ax.Toolbar    = [];
ax.SortMethod = 'childorder';   % skip depth-sorting -- irrelevant for 2D line data
disableDefaultInteractivity(ax);

% Static axes properties, set ONCE: channel labels/limits never change
% during navigation, so there is nothing to redo on every redraw.
ax.YTick      = tickPos;
ax.YTickLabel = tickLabels;
ax.YLim       = yLimits;
ax.XTickMode  = 'auto';

%%% --- Persistent graphics objects, created ONCE and updated via set() ---
% Dotted zero-line per channel: a fixed horizontal reference at each
% channel's baseline. XData spans a huge range so it always covers
% whatever the current view is, in any time unit -- it never needs
% touching again after this loop.
zeroLines = gobjects(nCh, 1);
for iCh = 1:nCh
    zeroLines(iCh) = line(ax, [-1e6 1e6], [offsets(iCh) offsets(iCh)], ...
        'Color', zeroLineGray, 'LineStyle', ':', 'LineWidth', 0.5, ...
        'HitTest', 'off', 'PickableParts', 'none');
end

% Each signal is ONE line spanning all channels (NaN gaps lift the pen
% between channels) -- created once, updated via XData/YData, never
% deleted/recreated. This is what a cla()+plot() approach was paying for
% every single redraw.
lineA = line(ax, NaN, NaN, 'Color', royalBlue, 'LineWidth', 0.5, ...
    'HitTest', 'off', 'PickableParts', 'none');
lineB = line(ax, NaN, NaN, 'Color', redColor, 'LineWidth', 0.5, ...
    'HitTest', 'off', 'PickableParts', 'none');

% Amplitude scale bar: a vertical reference line + label, updated (not
% recreated) every redraw since it depends on both ampScale and window.
scaleBarLine = line(ax, [NaN NaN], [NaN NaN], 'Color', 'k', 'LineWidth', 1.2, ...
    'HitTest', 'off', 'PickableParts', 'none');
scaleBarText = text(ax, NaN, NaN, '', 'VerticalAlignment', 'middle', ...
    'FontSize', 8, 'Color', 'k', 'HitTest', 'off');

xlabel(ax, 'Time (s)');

%%% --- Controls: single row, all sharing the same Y/height (aligned) ---
% Row lives at the very bottom, below the scrollbar (order bottom-to-top:
% buttons, scrollbar, hypnogram, EEG panel).
rowY = 0.02; rowH = 0.05/3; gap = 0.006;
xCur = 0.02;

pos = [xCur rowY 0.065 rowH]; xCur = xCur + pos(3) + gap;
uicontrol(fig, 'Style', 'pushbutton', 'String', char(9664) + " Prev", 'Units', 'normalized', ...
    'Position', pos, 'Callback', @(s,e) shiftEpoch(-1));

pos = [xCur rowY 0.065 rowH]; xCur = xCur + pos(3) + gap;
uicontrol(fig, 'Style', 'pushbutton', 'String', "Next " + char(9654), 'Units', 'normalized', ...
    'Position', pos, 'Callback', @(s,e) shiftEpoch(+1));

pos = [xCur rowY 0.045 rowH]; xCur = xCur + pos(3) + gap;
uicontrol(fig, 'Style', 'text', 'String', 'Epoch:', 'Units', 'normalized', ...
    'Position', pos, 'HorizontalAlignment', 'right', 'BackgroundColor', 'w');

pos = [xCur rowY 0.04 rowH]; xCur = xCur + pos(3) + gap;
editEpoch = uicontrol(fig, 'Style', 'edit', 'String', '1', 'Units', 'normalized', ...
    'Position', pos, 'Callback', @(s,e) jumpToEpoch(round(str2double(s.String))));
% NOTE: editEpoch is input-only -- it is deliberately NOT synced back to
% the current epoch on every redraw (that's what the title is for). That
% sync was one of the uicontrol touches that turned out to be expensive.

pos = [xCur rowY 0.065 rowH]; xCur = xCur + pos(3) + gap;
uicontrol(fig, 'Style', 'text', 'String', 'Window (s):', 'Units', 'normalized', ...
    'Position', pos, 'HorizontalAlignment', 'right', 'BackgroundColor', 'w');

pos = [xCur rowY 0.035 rowH]; xCur = xCur + pos(3) + gap;
uicontrol(fig, 'Style', 'edit', 'String', num2str(displaySec), 'Units', 'normalized', ...
    'Position', pos, 'Callback', @(s,e) setDisplaySec(str2double(s.String)));

pos = [xCur rowY 0.045 rowH]; xCur = xCur + pos(3) + gap;
uicontrol(fig, 'Style', 'text', 'String', 'Stage:', 'Units', 'normalized', ...
    'Position', pos, 'HorizontalAlignment', 'right', 'BackgroundColor', 'w');

pos = [xCur rowY 0.07 rowH]; xCur = xCur + pos(3) + gap;
popStage = uicontrol(fig, 'Style', 'popupmenu', 'String', stageNames, 'Units', 'normalized', ...
    'Position', pos);

pos = [xCur rowY 0.055 rowH]; xCur = xCur + pos(3) + gap;
uicontrol(fig, 'Style', 'pushbutton', 'String', 'First', 'Units', 'normalized', ...
    'Position', pos, 'Callback', @(s,e) jumpToStage('first'));

pos = [xCur rowY 0.08 rowH]; xCur = xCur + pos(3) + gap;
uicontrol(fig, 'Style', 'pushbutton', 'String', 'Next stage', 'Units', 'normalized', ...
    'Position', pos, 'Callback', @(s,e) jumpToStage('next'));

pos = [xCur rowY 0.075 rowH]; xCur = xCur + pos(3) + gap;
uicontrol(fig, 'Style', 'text', 'String', 'Time scale:', 'Units', 'normalized', ...
    'Position', pos, 'HorizontalAlignment', 'right', 'BackgroundColor', 'w');

pos = [xCur rowY 0.075 rowH]; xCur = xCur + pos(3) + gap;
popScale = uicontrol(fig, 'Style', 'popupmenu', 'String', {'seconds', 'minutes', 'hours'}, ...
    'Units', 'normalized', 'Position', pos, 'Callback', @(s,e) redraw());

pos = [xCur rowY 0.07 rowH]; xCur = xCur + pos(3) + gap;
btnDiff = uicontrol(fig, 'Style', 'togglebutton', 'String', 'Diff (d)', 'Units', 'normalized', ...
    'Position', pos, 'Callback', @(s,e) toggleDiff());

pos = [xCur rowY 0.04 rowH]; xCur = xCur + pos(3) + gap;
uicontrol(fig, 'Style', 'pushbutton', 'String', 'Amp -', 'Units', 'normalized', ...
    'Position', pos, 'Callback', @(s,e) scaleDown());

pos = [xCur rowY 0.04 rowH]; xCur = xCur + pos(3) + gap;
uicontrol(fig, 'Style', 'pushbutton', 'String', 'Amp +', 'Units', 'normalized', ...
    'Position', pos, 'Callback', @(s,e) scaleUp());

%%% --- Hypnogram strip (static, drawn once -- never touched in redraw) ---
% Shares panelX/panelW with the EEG panel and the scrollbar, so all three
% line up horizontally.
hypnoLabels = {'N3', 'N2', 'N1', 'R', 'W'};   % short form of stageNames, same order as stageValues
hypnoY = 0.12; hypnoH = 0.05;
axHypno = axes('Parent', fig, 'Units', 'normalized', ...
    'Position', [panelX, hypnoY, panelW, hypnoH]);
axHypno.YLim       = [min(stageValues)-0.5, max(stageValues)+0.5];
axHypno.YTick      = stageValues;
axHypno.YTickLabel = hypnoLabels;
axHypno.XTickMode  = 'auto';
axHypno.Box        = 'on';
axHypno.FontSize   = 7;
axHypno.TickLength = [0 0];
axHypno.Toolbar    = [];
disableDefaultInteractivity(axHypno);
xlabel(axHypno, 'Time (h)', 'FontSize', 7);

if ~isempty(scoringDigits)
    % XLim ends exactly at the last scored epoch -- not at totalDur, which
    % is the EEG recording's own length and isn't guaranteed to match the
    % scored range (that mismatch was leaving empty space at the end).
    hypnoEnd   = numel(scoringDigits) * scoreEpochSec;
    epochTimes = (0:numel(scoringDigits)-1) * scoreEpochSec;
    stairs(axHypno, [epochTimes, hypnoEnd] / 3600, double([scoringDigits(:); scoringDigits(end)]), ...
        'Color', [0.20 0.20 0.55], 'LineWidth', 1);
    axHypno.XLim = [0, hypnoEnd] / 3600;
else
    axHypno.XLim = [0, max(totalDur, displaySec)] / 3600;
    text(axHypno, mean(axHypno.XLim), mean(axHypno.YLim), 'No sleep scoring loaded', ...
        'HorizontalAlignment', 'center', 'FontSize', 8, 'Color', [0.5 0.5 0.5]);
end

%%% --- Custom scrollbar ---
% Drawn as patches inside a thin axes instead of a native uicontrol
% slider. On a GPU-less machine, moving a native Win32/Java slider turned
% out to cost 700-1300ms per redraw (confirmed by timing it in isolation
% from everything else), dwarfing the ~40ms the actual line drawing
% takes. A patch is rendered by the same painters pipeline as the traces
% and costs next to nothing by comparison.
sliderY = 0.05; sliderH = 0.035;
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

redraw();


%% ---- nested callbacks (share this function's workspace directly) ----

    function onKeyPress(~, evt)
        tCb = tic;
        switch evt.Key
            case 'leftarrow',  shiftEpoch(-1);
            case 'rightarrow', shiftEpoch(+1);
            case 'd',          toggleDiff();
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
        redraw();
        fprintf('[toggleDiff] TOTAL %6.1f ms\n\n', toc(tCb)*1000);
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
        tScaled = t / divisor;

        tStep = tic;
        Xall = repmat([tScaled, NaN], 1, nCh);
        if diffMode
            yMat = (EEG1data(:, idxStart:idxEnd) - EEG2data(:, idxStart:idxEnd)) * ampScale + offsets;
            Yall = reshape([yMat, nan(nCh, 1)].', 1, []);
            set(lineA, 'XData', Xall, 'YData', Yall, 'Color', diffColor, 'Visible', 'on');
            set(lineB, 'Visible', 'off');
        else
            yMat1 = EEG1data(:, idxStart:idxEnd) * ampScale + offsets;
            yMat2 = EEG2data(:, idxStart:idxEnd) * ampScale + offsets;
            Yall1 = reshape([yMat1, nan(nCh, 1)].', 1, []);
            Yall2 = reshape([yMat2, nan(nCh, 1)].', 1, []);
            set(lineA, 'XData', Xall, 'YData', Yall1, 'Color', royalBlue, 'Visible', 'on');
            set(lineB, 'XData', Xall, 'YData', Yall2, 'Color', redColor,  'Visible', 'on');
        end
        fprintf('  [lines]          %6.1f ms\n', toc(tStep)*1000);

        tStep = tic;
        target    = 0.3 * spacing;
        [~, bIdx] = min(abs(scaleCandidates * ampScale - target));
        barVal    = scaleCandidates(bIdx);
        barHeight = barVal * ampScale;
        barX      = tScaled(1) + 0.015 * (tScaled(end) - tScaled(1));
        barY0     = yLimits(2) - barHeight - spacing * 0.05;
        set(scaleBarLine, 'XData', [barX barX], 'YData', [barY0, barY0 + barHeight]);
        set(scaleBarText, 'Position', [barX + 0.01*(tScaled(end)-tScaled(1)), barY0 + barHeight/2, 0], ...
            'String', sprintf('%g %sV', barVal, char(956)));
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
        if diffMode, modeStr = '   |   DIFFERENCE (signal 1 - signal 2)'; end
        title(ax, sprintf('Epoch %d / %d   (%.1f-%.1f s)%s%s', ...
            curEpoch, nEpochs, t(1), t(end), stageStr, modeStr), 'Interpreter', 'none');
        fprintf('  [title]          %6.1f ms\n', toc(tStep)*1000);

        % Legend + diff button: only touched when diffMode actually
        % changes, not on every pure-navigation redraw.
        tStep = tic;
        if diffMode ~= lastDiffMode
            if diffMode
                legend(ax, lineA, {sprintf('Difference (%s - %s)', n1, n2)}, ...
                    'Location', 'northeast', 'Interpreter', 'none');
            else
                legend(ax, [lineA lineB], {n1, n2}, 'Location', 'northeast', 'Interpreter', 'none');
            end
            set(btnDiff, 'Value', diffMode);
            lastDiffMode = diffMode;
        end
        fprintf('  [legend/diffbtn] %6.1f ms\n', toc(tStep)*1000);

        tStep = tic;
        set(scrollThumb, 'XData', [windowStart, windowStart+displaySec, windowStart+displaySec, windowStart]);
        fprintf('  [scrollbar]      %6.1f ms\n', toc(tStep)*1000);

        tStep = tic;
        drawnow
        fprintf('  [drawnow]        %6.1f ms\n', toc(tStep)*1000);

        fprintf('  TOTAL redraw     %6.1f ms\n\n', toc(tTotal)*1000);
    end

end
