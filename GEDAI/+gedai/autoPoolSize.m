function n = autoPoolSize(opts)
% AUTOPOOLSIZE  Worker count for the FOOOF pool, scaled to the machine it runs on.
%   n = gedai.autoPoolSize()
%   n = gedai.autoPoolSize('MemoryPerWorkerGB', 2, 'Fraction', 0.25, 'MaxWorkers', 16)
%
%   Each worker is a separate MATLAB process with its own memory, and it competes with a
%   client that is already holding two full-night recordings plus the interpolated copies
%   run.eval_clean makes of them. A fixed worker count cannot be right across the range of
%   machines this pipeline runs on - a number chosen to leave headroom on a 64 GB box
%   throws away most of a 2 TB one, and a number chosen for the 2 TB box is what crowds
%   the small one. So budget a share of available memory instead and cap by core count,
%   which is what actually limits FOOOF once memory is not the binding constraint.
%
%   Fraction            share of currently-available memory the pool may occupy.
%                       Default 0.25 - the client needs the rest for the recordings.
%   MemoryPerWorkerGB   budget per worker. Default 2, comfortably above what a worker
%                       needs for its slice of the smoothed Fz spectra plus its own
%                       process overhead.
%   MaxWorkers          hard ceiling. Default: physical cores, further capped by the
%                       'Processes' cluster profile's own NumWorkers limit.
%
%   Worked examples, all with the defaults: a 64 GB / 16-core machine gets 7 workers, a
%   96 GB / 16-core machine gets 11, and a 2 TB machine is limited by its cores alone -
%   the memory term stops binding long before that, which is the point.
%
%   Returns at least 1. When memory cannot be measured (memory() is not available on every
%   platform and release) the memory term is dropped and only the core cap applies.

arguments
    opts.MemoryPerWorkerGB (1,1) double {mustBePositive} = 2
    opts.Fraction          (1,1) double {mustBePositive} = 0.25
    opts.MaxWorkers                                      = []
end

%%% Ceiling: cores, and whatever the cluster profile will actually hand out.
if isempty(opts.MaxWorkers)
    maxWorkers = feature('numcores');
else
    maxWorkers = opts.MaxWorkers;
end
try
    maxWorkers = min(maxWorkers, parcluster('Processes').NumWorkers);
catch
    %%% No Parallel Computing Toolbox, or no such profile - the core cap stands alone.
end

%%% Memory term. Skipped rather than guessed when memory() is unavailable, so a platform
%%% that cannot report leaves the core cap in charge instead of silently returning 1.
availableGB = NaN;
try
    m = memory;
    availableGB = m.MemAvailableAllArrays / 2^30;
catch
end

if isnan(availableGB)
    n = maxWorkers;
else
    n = floor((availableGB * opts.Fraction) / opts.MemoryPerWorkerGB);
    n = min(n, maxWorkers);
end
n = max(1, n);
end
