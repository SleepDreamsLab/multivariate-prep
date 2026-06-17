function ScoringDigits = scoreloader(filename)
% Handles multiple scoring types

[~, ~, ext] = fileparts(filename);

switch lower(ext)
    case '.csv'
        % Sleeptrip
        Scoring         = readmatrix(filename);
        ScoringDigits   = Scoring(:,1);  

    case '.json'
        % Scoringhero
        Scoring         = jsondecode(fileread(filename));
        ScoringDigits   = [Scoring{1}.digit];

    otherwise
        error('Unsupported file type: %s', ext);
end
end
