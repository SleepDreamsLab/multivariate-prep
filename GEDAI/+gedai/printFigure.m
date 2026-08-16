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
%   Returns true if a file was written.
%
% See also: print

if nargin < 3 || isempty(varargin), varargin = {'-r100'}; end

d = fileparts(fname);
if ~isempty(d) && ~exist(d, 'dir'), mkdir(d); end

try
    print(fig, fname, '-dpng', varargin{:});
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
    fprintf('gedai.printFigure: default renderer failed, saved with -painters instead.\n');
    ok = true;
    return
catch
end

warning('gedai:printFigure:failed', ...
    'Could not save %s (%s). Continuing without the figure.', fname, firstMsg);
ok = false;
end
