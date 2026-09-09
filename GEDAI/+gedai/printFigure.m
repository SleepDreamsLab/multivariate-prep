function ok = printFigure(fig, fname, varargin)
% PRINTFIGURE  print() a figure to PNG without ever throwing.
%
%   ok = gedai.printFigure(fig, fname)
%   ok = gedai.printFigure(fig, fname, '-r100')
%
%   Saving a figure must never be able to destroy a recording's worth of compute.
%   On a headless Linux box print() can fail with
%       "There was a problem while generating the output: Failed to export."
%   which, called from inside a pipeline stage's try block, throws away the Zapline
%   and CleanLine result that had already been computed but not yet written. Here the
%   failure costs a figure and nothing else.
%
%   The default renderer is tried first; on failure it retries with -painters, which
%   needs no OpenGL context and usually succeeds where the hardware path cannot get one.
%
%   A second failure mode never throws at all: on a machine with no GPU, print() goes
%   through SwiftShader (a CPU software OpenGL rasteriser - see the note in
%   run.eval_clean), and a figure heavy enough - many contour levels, per-axes
%   colormaps, alpha-blended patches, especially with other big figures still open and
%   competing for the same limited software context - can come back as a solid black
%   PNG with print() reporting success throughout. Every figure this draws has an
%   explicit white ('Color','w') background, so a written file with no near-white pixel
%   anywhere is never a real render; isBlankPng below catches that and routes it through
%   the same -painters retry as a thrown error would.
%
%   Returns true if a file was written.
%
% See also: print

if nargin < 3 || isempty(varargin), varargin = {'-r150'}; end

d = fileparts(fname);
if ~isempty(d) && ~exist(d, 'dir'), mkdir(d); end

try
    print(fig, fname, '-dpng', varargin{:});
    if isBlankPng(fname)
        error('gedai:printFigure:blank', 'default renderer produced a blank/black PNG');
    end
    ok = true;
    return
catch ME
    firstMsg = ME.message;
end

try
    %%% '-painters' rather than its modern spelling '-vector': the newer flag errors on
    %%% releases before R2022a, and this is the path that runs when things are already
    %%% going wrong. The Code Analyzer note about it is advisory.
    print(fig, fname, '-dpng', '-painters', varargin{:});
    if isBlankPng(fname)
        error('gedai:printFigure:blank', '-painters also produced a blank/black PNG');
    end
    fprintf('gedai.printFigure: default renderer failed, saved with -painters instead.\n');
    ok = true;
    return
catch
end

warning('gedai:printFigure:failed', ...
    'Could not save %s (%s). Continuing without the figure.', fname, firstMsg);
ok = false;
end

% -------------------------------------------------------------------------
function tf = isBlankPng(fname)
% Every figure gedai.printFigure is asked to save is drawn with a white figure
% background, so a genuine render always has pixels at or near 255 somewhere (the
% margins alone guarantee it). A file with nothing above a near-black ceiling is the
% SwiftShader failure this guards against, not a legitimately dark plot - none of ours
% are. Unreadable counts as blank too: either way the file is not fit to keep as-is.
try
    img = imread(fname);
    tf = max(img(:)) < 30;
catch
    tf = true;
end
end
