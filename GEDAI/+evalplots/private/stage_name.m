function name = stage_name(d)
switch d
    case  1,  name = 'Wake';
    case -1,  name = 'N1';
    case -2,  name = 'N2';
    case -3,  name = 'N3';
    case  0,  name = 'REM';
    case  5,  name = 'REM';
    otherwise, name = sprintf('Stage %d', d);
end
end
