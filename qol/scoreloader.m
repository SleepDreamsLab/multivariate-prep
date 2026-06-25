function scoringDigits = scoreloader(filename)
% Handles multiple scoring types

[~, ~, ext] = fileparts(filename);

switch lower(ext)
    case '.csv'
        % Sleeptrip
        Scoring         = readmatrix(filename);
        scoringDigits   = Scoring(:,1);  

        % Translate scoring
        scoremap = struct('From', [5  0  1  2  3 ], ...
                          'To',   [0  1  -1  -2  -3 ]);
        [tf, idx] = ismember(scoringDigits, scoremap.From);

        fprintf('\nOriginal scoring map\n')
        tabulate(scoringDigits)
        scoringDigits(tf) = scoremap.To(idx(tf));   
        
        fprintf('--> Relabeled scoring map\n')        
        tabulate(scoringDigits)

    case '.json'
        % Scoringhero
        Scoring         = jsondecode(fileread(filename));
        scoringDigits   = [Scoring{1}.digit];

    otherwise
        error('Unsupported file type: %s', ext);
end

% Make row vector
if ~isrow(scoringDigits)
    scoringDigits = scoringDigits';
end

end
