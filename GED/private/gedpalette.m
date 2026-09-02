function C = gedpalette(n)
% GEDPALETTE  Okabe-Ito: colours that survive colour blindness and greyscale.
%
%   Shared by plotged and plotgednight so a component keeps the same colour
%   across every panel of both figures. Yellow is last, being weakest on white.

base = [
    0.00 0.45 0.70      % blue
    0.84 0.37 0.00      % vermillion
    0.00 0.62 0.45      % bluish green
    0.80 0.47 0.65      % reddish purple
    0.34 0.71 0.91      % sky blue
    0.90 0.62 0.00      % orange
    0.35 0.35 0.35      % grey
    0.95 0.90 0.25];    % yellow
C = base(mod(0:n - 1, size(base, 1)) + 1, :);
end
