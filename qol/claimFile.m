function [guard, acquired, holder] = claimFile(lockFile, opts)
% CLAIMFILE  Take an exclusive, self-releasing claim on one unit of work.
%
%   [guard, acquired, holder] = claimFile(lockFile)
%   [...] = claimFile(lockFile, 'stalemin', 360, 'heartbeatmin', 5)
%
%   Lets several MATLAB sessions - on several machines - walk the same file list
%   without processing the same recording twice. The claim is a lock file created
%   with an atomic create-if-absent, so two sessions racing for the same recording
%   cannot both win: exactly one createNewFile call returns true. Checking isfile
%   first and then writing would not be enough; the gap between the two is exactly
%   where the race lives.
%
%   guard      onCleanup object. The claim is held for as long as it lives and is
%              released the moment it is cleared, overwritten, or its workspace is
%              destroyed - which covers the error and Ctrl-C paths, so a recording
%              that fails does not stay locked out forever. Keep it in a variable;
%              a discarded return value is destroyed at once and releases immediately.
%   acquired   false when someone else holds the claim (guard is then []).
%   holder     struct describing the current holder when acquired is false:
%              .host .pid .started .heartbeat (UTC ISO-8601 strings).
%
%   Name-value
%   ----------
%   stalemin       Minutes after which a claim whose heartbeat stopped is taken to be
%                  abandoned and is taken over. Must exceed the longest plausible
%                  runtime for one unit of work - see the heartbeat note below.
%                  Default: 360 (6 h).
%   heartbeatmin   Minutes between heartbeat writes. Default: 5.
%
%   Heartbeat note
%   --------------
%   The heartbeat runs on a MATLAB timer, and timer callbacks only fire when the
%   interpreter reaches its event queue - a long builtin or a parfor can hold that
%   off for minutes. Treat the heartbeat as a way to keep a slow job's claim fresh,
%   not as a liveness signal precise enough to justify a short stalemin. A machine
%   that crashes or reboots leaves its lock file behind: it is reclaimed only once
%   stalemin has passed, or when the lock file is deleted by hand.

arguments
    lockFile (1,:) char
    opts.stalemin     (1,1) double {mustBePositive} = 360
    opts.heartbeatmin (1,1) double {mustBePositive} = 5
end

guard = []; acquired = false; holder = emptyHolder();

lockDir = fileparts(lockFile);
if ~isempty(lockDir) && ~isfolder(lockDir), mkdir(lockDir); end

me = struct('host', localHost(), 'pid', feature('getpid'), ...
    'started', utcstamp(), 'heartbeat', utcstamp());

if ~createExclusive(lockFile)
    holder = readLock(lockFile);
    age    = claimAge(holder, lockFile);
    if age < opts.stalemin
        return                                  % someone else is on it
    end
    fprintf('[lock] claim on %s by %s (pid %d) is %.0f min stale - taking over\n', ...
        lockFile, holder.host, holder.pid, age);
    try delete(lockFile); catch, end
    if ~createExclusive(lockFile)               % another session beat us to the takeover
        holder = readLock(lockFile);
        return
    end
end

writeLock(lockFile, me);

%%% Heartbeat, so a job that outlives stalemin is not declared abandoned while it is
%%% still running. Failing to start one is not fatal - it only means the claim ages
%%% from its start time instead of from its last sign of life.
t = [];  %#ok<NASGU> replaced below when the timer starts
try
    period = max(1, opts.heartbeatmin * 60);
    t = timer('Name', 'gedai-claim-heartbeat', 'ExecutionMode', 'fixedSpacing', ...
        'Period', period, 'StartDelay', period, 'BusyMode', 'drop', ...
        'TimerFcn', @(~, ~) beat(lockFile, me));
    start(t);
catch ME
    warning('claimFile:noHeartbeat', ...
        'Could not start the claim heartbeat (%s); the claim ages from its start time.', ...
        ME.message);
    t = [];
end

guard    = onCleanup(@() release(lockFile, me, t));
acquired = true;
end

% -------------------------------------------------------------------------
function ok = createExclusive(f)
% Create f only if it does not exist yet, and report whether we were the creator.
% java.io.File.createNewFile does this in one atomic operation on both NTFS and SMB
% shares, which is what makes the claim race-free.
ok = false;
if usejava('jvm')
    try ok = logical(java.io.File(f).createNewFile()); catch, ok = false; end
    return
end
%%% JVM-less fallback (matlab -nojvm): check-then-create, so two sessions arriving in
%%% the same instant can both pass. Narrow enough to live with on a handful of
%%% machines, but the JVM path is the one that actually guarantees exclusivity.
if isfile(f), return; end
fid = fopen(f, 'w');
if fid >= 3, fclose(fid); ok = true; end
end

% -------------------------------------------------------------------------
function writeLock(lockFile, info)
fid = fopen(lockFile, 'w');
if fid < 3, return; end
fprintf(fid, '%s', jsonencode(info, 'PrettyPrint', true));
fclose(fid);
end

% -------------------------------------------------------------------------
function info = readLock(lockFile)
info = emptyHolder();
try
    got = jsondecode(fileread(lockFile));
    for f = intersect(fieldnames(info), fieldnames(got))'
        info.(f{1}) = got.(f{1});
    end
catch
    % Empty or half-written file: another session created it microseconds ago and has
    % not written its payload yet, or a write was interrupted. claimAge then falls back
    % to the file timestamp, which reads such a lock as fresh rather than abandoned.
end
end

% -------------------------------------------------------------------------
function mins = claimAge(holder, lockFile)
% Minutes since the holder last showed a sign of life.
if ~isempty(holder.heartbeat)
    try
        t = datetime(holder.heartbeat, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z''', ...
            'TimeZone', 'UTC');
        mins = minutes(datetime('now', 'TimeZone', 'UTC') - t);
        return
    catch
    end
end
%%% Fallback: file mtime. dir reports it in this machine's local time and `now` is
%%% local too, so the difference is right whatever timezone the other machine or the
%%% file server is set to.
mins = 0;
d = dir(lockFile);
if ~isempty(d)
    mins = minutes(datetime('now') - datetime(d(1).datenum, 'ConvertFrom', 'datenum'));
end
end

% -------------------------------------------------------------------------
function beat(lockFile, me)
try
    if ~isOurs(readLock(lockFile), me), return; end   % taken over - stop touching it
    me.heartbeat = utcstamp();
    writeLock(lockFile, me);
catch
end
end

% -------------------------------------------------------------------------
function release(lockFile, me, t)
if ~isempty(t) && isvalid(t)
    try stop(t); catch, end
    try delete(t); catch, end
end
%%% Only ever delete our own lock: if this claim was declared stale and taken over
%%% while we were still running, the file now belongs to another session.
if isOurs(readLock(lockFile), me)
    try
        delete(lockFile);
    catch ME
        warning('claimFile:releaseFailed', ...
            'Could not remove %s (%s) - delete it by hand, or that file stays claimed until it goes stale.', ...
            lockFile, ME.message);
    end
end
end

% -------------------------------------------------------------------------
function tf = isOurs(info, me)
tf = strcmp(info.host, me.host) && isequal(info.pid, me.pid);
end

% -------------------------------------------------------------------------
function info = emptyHolder()
info = struct('host', '', 'pid', NaN, 'started', '', 'heartbeat', '');
end

% -------------------------------------------------------------------------
function h = localHost()
h = getenv('COMPUTERNAME');
if isempty(h), h = getenv('HOSTNAME'); end
if isempty(h) && usejava('jvm')
    try h = char(java.net.InetAddress.getLocalHost().getHostName()); catch, end
end
if isempty(h), h = 'unknown-host'; end
end

% -------------------------------------------------------------------------
function s = utcstamp()
s = char(datetime('now', 'TimeZone', 'UTC', 'Format', 'yyyy-MM-dd''T''HH:mm:ss''Z'''));
end
