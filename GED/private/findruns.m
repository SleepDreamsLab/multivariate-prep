function runs = findruns(mask)
% FINDRUNS  Start/stop index pairs of each true run in a logical vector.
%
%   Returned as a 2 x n array so it can be walked with a for loop - each
%   iteration yields one column [start; stop], both indices inclusive.

mask   = mask(:)';
d      = diff([false mask false]);
starts = find(d == 1);
stops  = find(d == -1) - 1;
runs   = [starts; stops];
end
