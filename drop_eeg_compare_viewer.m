function drop_eeg_compare_viewer
%DROP_EEG_COMPARE_VIEWER
% Overlay two EEG recordings for direct visual comparison.
%
% Red  = signal 1
% Blue = signal 2
%
% Requires EEGLAB (pop_resample, pop_interp, pop_select, ...) already on
% the path.

%% Settings

WINSEC   = 30;
EPOCHSEC = 30;

file1       = '\\vs03.herseninstituut.knaw.nl\VS03-SandC-1\data\nin\data-drop\rawdata\derivatives\prep-ged\sub-drop0001\ses-t1\sub-drop0001_ses-t1_task-sleep_run-01_desc-zc2gedWakeBBAuto_eeg.set';   % first EEG file
file2       = '\\vs03.herseninstituut.knaw.nl\VS03-SandC-1\data\nin\data-drop\rawdata\derivatives\prep-ged\sub-drop0001\ses-t1\sub-drop0001_ses-t1_task-sleep_run-01_desc-zc2gedWakeBBAutoPlusFSAutoPlus_eeg.set';   % second EEG file
scoringFile = '\\vs03.herseninstituut.knaw.nl\VS03-SandC-1\data\nin\data-drop\derivatives\scores\final\sub-drop0001_ses-t1_task-sleep_run-01_eeg.csv';   % sleep scoring file, '' if none

net = 'EGI256';   % net type, passed to chans1020

default1020Labels = [
    "Fp1"
    "Fpz"
    "Fp2"
    "F7"
    "F3"
    "Fz"
    "F4"
    "F8"
    "T7"
    "C3"
    "Cz"
    "C4"
    "T8"
    "P7"
    "P3"
    "Pz"
    "P4"
    "P8"
    "O1"
    "Oz"
    "O2"
    ];

%% Load + prepare EEG
% Same import path as EEGCompareViewer.m: fast_eeg_import, then interpolate
% missing channels and keep only the 10-20 electrodes via chans1020.

EEG1 = fast_eeg_import(file1);
EEG2 = fast_eeg_import(file2);

EEG1 = prepEEG1020(EEG1, net);
EEG2 = prepEEG1020(EEG2, net);

%%% --- Align channels between the two files (label match, stable order) ---
labels1 = {EEG1.chanlocs.labels};
labels2 = {EEG2.chanlocs.labels};
[~, ia, ib] = intersect(labels1, labels2, 'stable');
if numel(ia) < numel(labels1) || numel(ib) < numel(labels2)
    warning('drop_eeg_compare_viewer:chanMismatch', ...
        'Channel labels differ between files; using %d common channels.', numel(ia));
end
EEG1.data     = EEG1.data(ia, :);
EEG1.chanlocs = EEG1.chanlocs(ia);
EEG2.data     = EEG2.data(ib, :);

if EEG1.srate ~= EEG2.srate
    warning('drop_eeg_compare_viewer:srateMismatch', ...
        'Sampling rates differ (%.3f vs %.3f Hz); using file 1''s rate.', EEG1.srate, EEG2.srate);
end

data1 = EEG1.data;
data2 = EEG2.data;
srate = EEG1.srate;
pnts  = min(size(data1, 2), size(data2, 2));
data1 = data1(:, 1:pnts);
data2 = data2(:, 1:pnts);

channelLabels = {EEG1.chanlocs.labels};

%% Select default 10-20 channels to display

[selectedChannels, matched1020, missing1020] = ...
    findRequestedChannels(channelLabels, default1020Labels);

if isempty(selectedChannels)
    selectedChannels = 1:numel(channelLabels);
    matched1020 = string(channelLabels);
    missing1020 = strings(0, 1);
end

defaultChannels = selectedChannels;

fprintf('\nDefault display channels:\n');
fprintf('%s\n', strjoin(matched1020, ', '));

if ~isempty(missing1020)
    fprintf('\nRequested 10-20 channels not found:\n');
    fprintf('%s\n', strjoin(missing1020, ', '));
end

%% Load sleep scoring

stageNames  = {'N3', 'N2', 'N1', 'REM', 'Wake'};
stageValues = [-3, -2, -1, 0, 1];

stageCodes = [];
if ~isempty(scoringFile)
    stageCodes = scoreloader(scoringFile);
    fprintf('\nSleep scoring loaded from:\n%s\n', scoringFile);
else
    fprintf('\nNo sleep scoring loaded. Stage navigation is disabled.\n');
end

%% Viewer state

epochSamp   = round(EPOCHSEC * srate);
nDataEpochs = max(1, floor(pnts / epochSamp));
winSamp     = min(round(WINSEC * srate), pnts);
maxStart    = max(1, pnts - winSamp + 1);
t0          = 1;

%% Create figure

[~, n1] = fileparts(file1);
[~, n2] = fileparts(file2);

fig = figure( ...
    'Name', sprintf('%s  vs  %s', n1, n2), ...
    'NumberTitle', 'off', ...
    'Color', 'w', ...
    'Position', [70 60 1450 850], ...
    'MenuBar', 'none', ...
    'WindowScrollWheelFcn', @onMouseWheel, ...
    'WindowKeyPressFcn', @onKeyPress);

ax = axes( ...
    'Parent', fig, ...
    'Position', [0.055 0.10 0.92 0.76]);

uicontrol(fig, ...
    'Style', 'pushbutton', ...
    'Units', 'normalized', ...
    'Position', [0.055 0.94 0.105 0.035], ...
    'String', 'Select channels', ...
    'FontSize', 9, ...
    'Callback', @chooseChannels);

uicontrol(fig, ...
    'Style', 'pushbutton', ...
    'Units', 'normalized', ...
    'Position', [0.168 0.94 0.105 0.035], ...
    'String', '10-20 defaults', ...
    'FontSize', 9, ...
    'Callback', @restoreChannelDefaults);

statusText = uicontrol(fig, ...
    'Style', 'text', ...
    'Units', 'normalized', ...
    'Position', [0.29 0.935 0.68 0.045], ...
    'BackgroundColor', 'w', ...
    'HorizontalAlignment', 'left', ...
    'FontSize', 9, ...
    'String', '');

uicontrol(fig, ...
    'Style', 'text', ...
    'Units', 'normalized', ...
    'Position', [0.055 0.892 0.045 0.030], ...
    'BackgroundColor', 'w', ...
    'HorizontalAlignment', 'left', ...
    'String', 'Epoch:');

epochEdit = uicontrol(fig, ...
    'Style', 'edit', ...
    'Units', 'normalized', ...
    'Position', [0.100 0.895 0.055 0.032], ...
    'BackgroundColor', 'w', ...
    'String', '1', ...
    'Callback', @jumpToEpoch);

uicontrol(fig, ...
    'Style', 'pushbutton', ...
    'Units', 'normalized', ...
    'Position', [0.162 0.893 0.070 0.035], ...
    'String', 'Jump', ...
    'Callback', @jumpToEpoch);

uicontrol(fig, ...
    'Style', 'text', ...
    'Units', 'normalized', ...
    'Position', [0.255 0.892 0.040 0.030], ...
    'BackgroundColor', 'w', ...
    'HorizontalAlignment', 'left', ...
    'String', 'Stage:');

stagePopup = uicontrol(fig, ...
    'Style', 'popupmenu', ...
    'Units', 'normalized', ...
    'Position', [0.297 0.893 0.075 0.035], ...
    'BackgroundColor', 'w', ...
    'String', stageNames, ...
    'Value', 3);

firstStageButton = uicontrol(fig, ...
    'Style', 'pushbutton', ...
    'Units', 'normalized', ...
    'Position', [0.380 0.893 0.060 0.035], ...
    'String', 'First', ...
    'Callback', @(~,~) jumpToStage('first'));

previousStageButton = uicontrol(fig, ...
    'Style', 'pushbutton', ...
    'Units', 'normalized', ...
    'Position', [0.447 0.893 0.075 0.035], ...
    'String', 'Previous', ...
    'Callback', @(~,~) jumpToStage('previous'));

nextStageButton = uicontrol(fig, ...
    'Style', 'pushbutton', ...
    'Units', 'normalized', ...
    'Position', [0.529 0.893 0.060 0.035], ...
    'String', 'Next', ...
    'Callback', @(~,~) jumpToStage('next'));

scoreText = uicontrol(fig, ...
    'Style', 'text', ...
    'Units', 'normalized', ...
    'Position', [0.610 0.888 0.36 0.040], ...
    'BackgroundColor', 'w', ...
    'HorizontalAlignment', 'left', ...
    'FontSize', 8, ...
    'String', '');

% NOTE: this is a native uicontrol slider (unlike EEGCompareViewer.m,
% which replaced its slider with a hand-drawn patch scrollbar after
% finding the native control's repaint cost 700-1300ms per redraw on a
% GPU-less machine). Kept here to match the requested GUI layout as-is.
timeSlider = uicontrol(fig, ...
    'Style', 'slider', ...
    'Units', 'normalized', ...
    'Position', [0.055 0.025 0.65 0.025], ...
    'Min', 1, ...
    'Max', max(2, maxStart), ...
    'Value', 1, ...
    'Callback', @setTimeFromSlider);

uicontrol(fig, ...
    'Style', 'pushbutton', ...
    'Units', 'normalized', ...
    'Position', [0.72 0.017 0.07 0.042], ...
    'String', 'Previous', ...
    'Callback', @(~,~) shiftTime(-winSamp));

uicontrol(fig, ...
    'Style', 'pushbutton', ...
    'Units', 'normalized', ...
    'Position', [0.80 0.017 0.07 0.042], ...
    'String', 'Next', ...
    'Callback', @(~,~) shiftTime(winSamp));

timeText = uicontrol(fig, ...
    'Style', 'text', ...
    'Units', 'normalized', ...
    'Position', [0.88 0.017 0.10 0.042], ...
    'BackgroundColor', 'w', ...
    'HorizontalAlignment', 'left', ...
    'FontSize', 9, ...
    'String', '');

if maxStart > 1
    sliderRange = maxStart - 1;
    smallStep = min(max(round(srate), 1) / sliderRange, 1);
    largeStep = min(winSamp / sliderRange, 1);
    timeSlider.SliderStep = [smallStep largeStep];
else
    set(timeSlider, 'Enable', 'off');
end

if isempty(stageCodes)
    set([stagePopup firstStageButton previousStageButton nextStageButton], ...
        'Enable', 'off');
    set(scoreText, 'String', 'No sleep scoring loaded');
else
    [~, scoreName, scoreExt] = fileparts(scoringFile);
    set(scoreText, 'String', sprintf('Scoring: %s%s', scoreName, scoreExt));
end

updateDisplay();

%% Callbacks

    function chooseChannels(~, ~)
        [selection, confirmed] = listdlg( ...
            'ListString', channelLabels, ...
            'SelectionMode', 'multiple', ...
            'InitialValue', selectedChannels, ...
            'Name', 'Select EEG channels', ...
            'PromptString', 'Channels to display:', ...
            'ListSize', [270 520]);

        if ~confirmed || isempty(selection)
            return
        end

        selectedChannels = selection(:).';
        updateDisplay();
    end

    function restoreChannelDefaults(~, ~)
        selectedChannels = defaultChannels;
        updateDisplay();
    end

    function jumpToEpoch(~, ~)
        requestedEpoch = str2double(get(epochEdit, 'String'));

        if ~isfinite(requestedEpoch) || ...
                requestedEpoch < 1 || ...
                requestedEpoch > nDataEpochs || ...
                requestedEpoch ~= round(requestedEpoch)

            warndlg(sprintf( ...
                'Enter a whole-number epoch between 1 and %d.', ...
                nDataEpochs), ...
                'Invalid epoch');

            currentEpoch = floor((t0 - 1) / epochSamp) + 1;
            set(epochEdit, 'String', num2str(currentEpoch));
            return
        end

        targetSample = (requestedEpoch - 1) * epochSamp + 1;
        setTime(targetSample);
    end

    function jumpToStage(direction)
        if isempty(stageCodes)
            return
        end

        targetCode = stageValues(get(stagePopup, 'Value'));
        matchingEpochs = find(stageCodes == targetCode);
        matchingEpochs = matchingEpochs( ...
            matchingEpochs >= 1 & matchingEpochs <= nDataEpochs);

        if isempty(matchingEpochs)
            warndlg(sprintf( ...
                'No %s epochs were found in this recording.', ...
                stageNames{get(stagePopup, 'Value')}), ...
                'Stage not found');
            return
        end

        currentEpoch = floor((t0 - 1) / epochSamp) + 1;

        switch direction
            case 'first'
                targetEpoch = matchingEpochs(1);

            case 'previous'
                candidate = matchingEpochs(matchingEpochs < currentEpoch);
                if isempty(candidate)
                    targetEpoch = matchingEpochs(end);
                else
                    targetEpoch = candidate(end);
                end

            case 'next'
                candidate = matchingEpochs(matchingEpochs > currentEpoch);
                if isempty(candidate)
                    targetEpoch = matchingEpochs(1);
                else
                    targetEpoch = candidate(1);
                end
        end

        setTime((targetEpoch - 1) * epochSamp + 1);
    end

    function setTimeFromSlider(source, ~)
        setTime(round(source.Value));
    end

    function shiftTime(numberOfSamples)
        setTime(t0 + numberOfSamples);
    end

    function setTime(newStart)
        t0 = min(max(1, round(newStart)), maxStart);

        if maxStart > 1
            set(timeSlider, 'Value', t0);
        end

        updateDisplay();
    end

    function onMouseWheel(~, event)
        step = round(5 * srate);
        shiftTime(event.VerticalScrollCount * step);
    end

    function onKeyPress(~, event)
        switch event.Key
            case 'rightarrow'
                shiftTime(winSamp);

            case 'leftarrow'
                shiftTime(-winSamp);
        end
    end

    function updateDisplay()
        i0 = t0;
        i1 = min(pnts, t0 + winSamp - 1);

        sampleIndices = i0:i1;
        time = (sampleIndices - 1) / srate;

        currentEpoch = floor((i0 - 1) / epochSamp) + 1;
        currentStage = stageNameForEpoch(stageCodes, currentEpoch, stageNames, stageValues);

        set(epochEdit, 'String', num2str(currentEpoch));

        signal1Shown = double(data1(selectedChannels, sampleIndices));
        signal2Shown = double(data2(selectedChannels, sampleIndices));

        combinedData = [signal1Shown signal2Shown];
        channelSD = std(combinedData, 0, 2, 'omitnan');
        channelSD = channelSD(isfinite(channelSD) & channelSD > 0);

        if isempty(channelSD)
            displayScale = 1;
        else
            displayScale = median(channelSD);
        end

        spacing = 6 * displayScale;

        plotOverlayEEG( ...
            ax, time, signal1Shown, signal2Shown, ...
            channelLabels(selectedChannels), spacing);

        title(ax, sprintf( ...
            ['%s (red) and %s (blue)' ...
             ' | epoch %d | %s | %.3f-%.3f h'], ...
            n1, n2, currentEpoch, currentStage, ...
            time(1)/3600, time(end)/3600), ...
            'Interpreter', 'none', ...
            'FontWeight', 'bold');

        xlabel(ax, 'Time (s)');

        set(statusText, 'String', sprintf( ...
            'Epoch %d/%d: %s   |   %d channel(s)', ...
            currentEpoch, nDataEpochs, currentStage, numel(selectedChannels)));

        set(timeText, 'String', sprintf( ...
            '%.2f-%.2f h', time(1)/3600, time(end)/3600));

        fig.Name = sprintf( ...
            '%s vs %s | %.2f-%.2f h', n1, n2, time(1)/3600, time(end)/3600);

        drawnow limitrate
    end

end


%% ===================== LOCAL FUNCTIONS =====================

function EEG = prepEEG1020(EEG, net)
    %%% Interpolate missing channels, then keep only 10-20 electrodes
    EEG = pop_resample(EEG, 100);
    try; EEG = pop_interp(EEG, EEG.urchanlocs, 'spherical'); end
    EEG = chans1020(EEG, false, 'add_eog', true, 'net', net);
end


function plotOverlayEEG(ax, time, signal1Data, signal2Data, labels, spacing)

cla(ax);

nChannels = size(signal1Data, 1);

if ~isfinite(spacing) || spacing <= 0
    spacing = 1;
end

plotOffsets = (nChannels-1:-1:0)' * spacing;

signal1Plot = signal1Data + plotOffsets;
signal2Plot = signal2Data + plotOffsets;

hold(ax, 'on');

% Signal 1 is drawn first, so it remains behind
plot(ax, time, signal1Plot.', ...
    'Color', [0.85 0.15 0.15], ...
    'LineWidth', 0.45);

% Signal 2 is drawn on top
plot(ax, time, signal2Plot.', ...
    'Color', [0.00 0.25 0.75], ...
    'LineWidth', 0.65);

hold(ax, 'off');

tickPositions = flipud(plotOffsets);
tickLabelsOut = flipud(labels(:));

set(ax, ...
    'YTick', tickPositions, ...
    'YTickLabel', tickLabelsOut, ...
    'TickLabelInterpreter', 'none', ...
    'FontSize', 8);

if numel(time) > 1
    xlim(ax, [time(1) time(end)]);
end

ylim(ax, [-spacing plotOffsets(1)+spacing]);

ylabel(ax, 'Channel');
box(ax, 'on');
grid(ax, 'on');

end


function [indices, matchedLabels, missingLabels] = ...
    findRequestedChannels(channelLabels, requestedLabels)

available = string(channelLabels(:));
requested = string(requestedLabels(:));

indices       = [];
matchedLabels = strings(0, 1);
missingLabels = strings(0, 1);

for k = 1:numel(requested)

    idx = find(strcmpi(strtrim(available), ...
        strtrim(requested(k))), 1);

    if isempty(idx)

        escapedLabel = regexptranslate( ...
            'escape', char(requested(k)));

        pattern = [ ...
            '(^|[^A-Za-z0-9])' ...
            escapedLabel ...
            '([^A-Za-z0-9]|$)'];

        matches = regexpi( ...
            cellstr(available), pattern, 'once');

        idx = find(~cellfun(@isempty, matches), 1);
    end

    if isempty(idx)
        missingLabels(end+1, 1) = requested(k);
    elseif ~ismember(idx, indices)
        indices(end+1) = idx;
        matchedLabels(end+1, 1) = requested(k);
    end
end

end


function stageName = stageNameForEpoch(stageCodes, epochNumber, stageNames, stageValues)

if isempty(stageCodes) || ...
        epochNumber < 1 || ...
        epochNumber > numel(stageCodes)

    stageName = 'no score';
    return
end

code = stageCodes(epochNumber);

if ~isfinite(code)
    stageName = 'unscored';
    return
end

hit = find(stageValues == code, 1);
if isempty(hit)
    stageName = sprintf('stage code %g', code);
else
    stageName = stageNames{hit};
end

end
