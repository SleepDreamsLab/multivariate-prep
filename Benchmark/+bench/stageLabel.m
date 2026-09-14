function label = stageLabel(digit)
% STAGELABEL  Pipeline stage digit (after scoreloader) to a readable label.
switch digit
    case -3
        label = 'N3';
    case -2
        label = 'N2';
    case -1
        label = 'N1';
    case 0
        label = 'REM';
    case 1
        label = 'Wake';
    otherwise
        label = sprintf('stage%d', digit);
end
end
