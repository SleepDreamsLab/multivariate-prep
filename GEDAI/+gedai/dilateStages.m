function [digits, nChanged] = dilateStages(digits, opts)
% DILATESTAGES  Grow sleep stages into their neighbouring epochs, for cleaning only.
%
%   digits = gedai.dilateStages(digits, 'from', [-2 -3 0], 'into', [1 -1], 'n', 1)
%   [digits, nChanged] = gedai.dilateStages(...)
%
%   Why: scoring assigns one stage to a whole 30-s epoch, so an epoch scored Wake or N1
%   at a transition can still contain genuine sleep slow waves. GEDAI cleans each stage
%   group with its own strength, and a stronger wake setting will treat those slow waves
%   as artefact. Relabelling the boundary epochs as sleep makes them get the gentler
%   sleep cleaning instead.
%
%   This is a CLEANING-TIME relabelling. Keep the original scoring for analysis and for
%   the evaluation figures - only the vector handed to GEDAI should be dilated, or your
%   stage contrasts will silently include relabelled epochs.
%
%   Run it BEFORE gedai.killN1: with N1 in 'into', an N1 epoch touching a sleep epoch is
%   claimed by that sleep stage first, and killN1 then only splits whatever N1 is left.
%
%   Name-value
%   ----------
%   from       Stages that grow. Order sets priority: when an epoch borders two of them,
%              the one listed first claims it.                 (default [-2 -3 0])
%   into       Stages that may be overwritten. Everything else is left alone - notably,
%              do not let one analysed stage overwrite another (default [1 -1])
%   n          How many epochs to grow on each side            (default 1)
%   direction  'both'      overwrite the epochs before and after a 'from' run
%              'forward'   only the epoch after  (later in time)
%              'backward'  only the epoch before (earlier in time)  (default 'forward')
%
%   Outputs
%   -------
%   digits     Relabelled scoring, same orientation as the input.
%   nChanged   Number of epochs relabelled - worth recording in the sidecar, since these
%              epochs also join the covariance GEDAI estimates its spatial filters from.
%
%   Example (defaults: forward only)
%   -------
%   [1 1 -1 -2 -2 1]  ->  [1 1 -1 -2 -2 -2]     the epoch after the N2 run is claimed;
%                                               the N1 before it is not
%   with direction='both':  ->  [1 1 -2 -2 -2 -2]
%
% See also: gedai.killN1

arguments
    digits            {mustBeVector, mustBeNumeric}
    opts.from   (1,:) double = [-2 -3 0]
    opts.into   (1,:) double = [1 -1]
    opts.n      (1,1) double {mustBeInteger, mustBeNonnegative} = 1
    opts.direction {mustBeMember(opts.direction, {'both', 'forward', 'backward'})} = 'forward'
end

wasColumn = iscolumn(digits);
orig      = digits(:)';
digits    = orig;

%%% Epochs claimed by this call are locked, so a later 'from' stage (or a later step)
%%% cannot take them again. Without this, priority would depend on iteration order.
locked = false(size(digits));

doFwd = any(strcmp(opts.direction, {'both', 'forward'}));
doBwd = any(strcmp(opts.direction, {'both', 'backward'}));

for iStep = 1:opts.n
    %%% Neighbours are read from the state at the start of the step, so each step grows
    %%% every stage by exactly one epoch rather than letting the first one run away.
    prev      = digits;
    claimable = ismember(prev, opts.into) & ~locked;

    for s = opts.from
        isS  = (prev == s);
        grow = false(size(prev));
        if doFwd, grow = grow | [false, isS(1:end-1)]; end   % epoch after  an s epoch
        if doBwd, grow = grow | [isS(2:end), false];  end    % epoch before an s epoch

        take            = grow & claimable;
        digits(take)    = s;
        locked(take)    = true;
        claimable(take) = false;      % first stage in 'from' wins the tie
    end
end

nChanged = nnz(digits ~= orig);
if nChanged > 0
    fprintf('gedai.dilateStages: relabelled %d epoch(s) (%s, n=%d).\n', ...
        nChanged, opts.direction, opts.n);
end

if wasColumn, digits = digits(:); end
end
