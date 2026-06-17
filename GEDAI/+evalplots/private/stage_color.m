function col = stage_color(d)
switch d
    case  1,  col = [0.82 0.10 0.15];   % crimson  – Wake
    case -1,  col = [0.20 0.55 0.80];   % sky blue – N1
    case -2,  col = [0.20 0.65 0.30];   % green    – N2
    case -3,  col = [0.50 0.20 0.75];   % purple   – N3
    case  0,  col = [0.05 0.37 0.73];   % blue     – REM
    case  5,  col = [0.05 0.37 0.73];   % blue     – REM
    otherwise, col = [0.40 0.40 0.40];
end
end
