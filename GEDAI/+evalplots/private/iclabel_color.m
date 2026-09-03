function [rgb, names] = iclabel_color(classIdx)
% ICLABEL_COLOR  Colour of an ICLabel class, in ICLabel's own class order.
%
%   rgb          = iclabel_color(classIdx)
%   [rgb, names] = iclabel_color()
%
%   classIdx indexes the ICLabel class list - 1=Brain, 2=Muscle, 3=Eye, 4=Heart,
%   5=Line Noise, 6=Channel Noise, 7=Other - and may be a vector, in which case
%   rgb comes back one row per element. An index outside the list (or an empty
%   one, i.e. a component with no classification) gets the neutral grey used for
%   "unknown", so a caller never has to special-case an unlabelled decomposition.
%
%   Same palette as ica.plot.amica_iclabel_bars, deliberately: a component that
%   is orange in the ICLabel bar chart should be orange here too, or reading the
%   two side by side means re-learning the code.

persistent tbl
if isempty(tbl)
    tbl = [
        0.20 0.55 0.30;   % Brain
        0.85 0.55 0.10;   % Muscle
        0.75 0.20 0.65;   % Eye
        0.80 0.20 0.20;   % Heart
        0.30 0.55 0.85;   % Line Noise
        0.55 0.40 0.25;   % Channel Noise
        0.55 0.55 0.55];  % Other
end

names = {'Brain', 'Muscle', 'Eye', 'Heart', 'Line Noise', 'Channel Noise', 'Other'};

unknown = [0.35 0.35 0.35];
if nargin < 1 || isempty(classIdx)
    rgb = unknown;
    return
end

classIdx = classIdx(:);
rgb      = repmat(unknown, numel(classIdx), 1);
ok       = isfinite(classIdx) & classIdx >= 1 & classIdx <= size(tbl, 1);
rgb(ok, :) = tbl(classIdx(ok), :);
end
