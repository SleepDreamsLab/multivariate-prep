function ScoringDigits = scoreloader(filename)
% Handles multiple scoring types

[~, ~, ext] = fileparts(filename);

switch lower(ext)
    case '.csv'
        % Sleeptrip
        Scoring         = readmatrix(filename);
        ScoringDigits   = Scoring(:,1);  

        % Translate scoring
        scoremap = struct('From', [5  0  -1  -2  -3 ], ...
                          'To',   [0  1  -1  -2  -3 ]);
        [tf, idx] = ismember(ScoringDigits, opts.ScoringMap.From);
        ScoringDigits(tf) = opts.ScoringMap.To(idx(tf));         

    case '.json'
        % Scoringhero
        Scoring         = jsondecode(fileread(filename));
        ScoringDigits   = [Scoring{1}.digit];

    otherwise
        error('Unsupported file type: %s', ext);
end
end
