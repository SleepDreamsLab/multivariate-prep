function scoringDigits = scoreloader(filename)
% Handles multiple scoring types

[~, ~, ext] = fileparts(filename);

scoremap = struct('From', [5  0  1  2  3], ...
                  'To',   [0  1  -1  -2  -3]);

switch lower(ext)
    case '.csv'
        % Sleeptrip
        Scoring         = readmatrix(filename);
        scoringDigits   = Scoring(:,1);

        fprintf('\nOriginal scoring map\n')
        tabulate(scoringDigits)
        [tf, idx] = ismember(scoringDigits, scoremap.From);
        scoringDigits(tf) = scoremap.To(idx(tf));
        fprintf('--> Relabeled scoring map\n')
        tabulate(scoringDigits)

    case '.json'
        % Scoringhero
        Scoring       = jsondecode(fileread(filename));
        scoringDigits = [Scoring{1}.digit];

    case '.xml'
        % Compumedics / Nox EventExport XML
        [~, hypnogram, ~] = extract_sleep_stages(filename);
        scoringDigits = hypnogram;

        fprintf('\nOriginal scoring map\n')
        tabulate(scoringDigits(~isnan(scoringDigits)))
        [tf, idx] = ismember(scoringDigits, scoremap.From);
        scoringDigits(tf) = scoremap.To(idx(tf));
        fprintf('--> Relabeled scoring map\n')
        tabulate(scoringDigits(~isnan(scoringDigits)))

    otherwise
        error('Unsupported file type: %s', ext);
end

% Make row vector
if ~isrow(scoringDigits)
    scoringDigits = scoringDigits';
end

end







% =========================================================================
function [stageTable, hypnogram, epochSec] = extract_sleep_stages(xmlFile)
% EXTRACT_SLEEP_STAGES  Parse sleep-stage events from a Compumedics/Nox-style
% "EventExport" XML (<Events><Event><Type>SLEEP-...</Type>...) and return
% them as a table plus a per-epoch hypnogram vector.

if nargin < 1 || ~isfile(xmlFile)
    error('extract_sleep_stages:badInput', 'xmlFile not found: %s', xmlFile);
end

doc        = xmlread(xmlFile);
eventNodes = doc.getElementsByTagName('Event');
nEvents    = eventNodes.getLength();

stageStr = strings(nEvents, 1);
startStr = strings(nEvents, 1);
stopStr  = strings(nEvents, 1);
keep     = false(nEvents, 1);

for i = 1:nEvents
    ev   = eventNodes.item(i - 1);
    type = getChildText(ev, 'Type');
    if startsWith(type, "SLEEP-")
        stageStr(i) = type;
        startStr(i) = getChildText(ev, 'StartTime');
        stopStr(i)  = getChildText(ev, 'StopTime');
        keep(i)     = true;
    end
end

if ~any(keep)
    error('extract_sleep_stages:noEvents', 'No SLEEP-* events found in %s', xmlFile);
end

stageStr = stageStr(keep);
startStr = startStr(keep);
stopStr  = stopStr(keep);

fmt       = 'yyyy-MM-dd''T''HH:mm:ss.SSSSSS';
StartTime = datetime(startStr, 'InputFormat', fmt);
StopTime  = datetime(stopStr,  'InputFormat', fmt);

stageTable             = table(StartTime, StopTime, stageStr, 'VariableNames', {'StartTime', 'StopTime', 'Stage'});
stageTable.DurationSec = seconds(stageTable.StopTime - stageTable.StartTime);
stageTable.StageCode   = stageCode(stageTable.Stage);
stageTable             = sortrows(stageTable, 'StartTime');
stageTable             = stageTable(:, {'StartTime', 'StopTime', 'DurationSec', 'Stage', 'StageCode'});

epochSec = mode(stageTable.DurationSec);
t0       = stageTable.StartTime(1);
epochIdx = round(seconds(stageTable.StartTime - t0) / epochSec) + 1;

hypnogram           = nan(max(epochIdx), 1);
hypnogram(epochIdx) = stageTable.StageCode;

end

% -------------------------------------------------------------------------
function txt = getChildText(parentNode, tagName)
txt      = "";
children = parentNode.getChildNodes();
for k = 0:children.getLength() - 1
    node = children.item(k);
    if node.getNodeType() == node.ELEMENT_NODE && strcmp(char(node.getNodeName()), tagName)
        firstChild = node.getFirstChild();
        if ~isempty(firstChild)
            txt = string(char(firstChild.getData()));
        end
        return
    end
end
end

% -------------------------------------------------------------------------
function code = stageCode(stageStr)
code = nan(size(stageStr));
code(stageStr == "SLEEP-S0")  = 0;
code(stageStr == "SLEEP-S1")  = 1;
code(stageStr == "SLEEP-S2")  = 2;
code(stageStr == "SLEEP-S3")  = 3;
code(stageStr == "SLEEP-S4")  = 3;
code(stageStr == "SLEEP-REM") = 5;
% SLEEP-UNSCORED (and anything unrecognized) stays NaN
end
