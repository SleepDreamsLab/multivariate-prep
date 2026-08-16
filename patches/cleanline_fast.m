function EEG = cleanline_fast(EEG, varargin)
% Faster drop-in for the legacy pop_cleanline() channel loop.
%
% Usage:
%   EEG = cleanline_fast(EEG, 'linefreqs', 50, 'winsize', 4, 'winstep', 2);
%
% Two speedups over the per-channel loop in cleanline():
%   1. Channels are distributed over a parpool with parfor. This is the bulk
%      of the gain.
%   2. rmlinesmovingwinc() is already vectorised over columns, so channels are
%      handed to it in small blocks (one taper set and one batched fft per
%      window instead of one per channel). Worth ~1.3x on its own; blocks much
%      larger than a handful of channels get slower again, hence blocksize 4.
% The spectral power computation (cleanline's 'computepower', two extra
% mtspectrumsegc() passes per channel) is dropped, since it only feeds the
% Sorig/Sclean return values.
%
% Options (name/value, defaults match pop_cleanline):
%   chanlist    channels to clean                        (all)
%   linefreqs   line frequencies to remove               (60)
%   p           p-value for sinusoid detection           (0.01)
%   bandwidth   spectral bandwidth of a line (Hz)        (2)
%   winsize     sliding window length (sec)              (4)
%   winstep     sliding window step (sec)                (1)
%   tau         overlap smoothing factor                 (100)
%   blocksize   channels per parfor block                (4)
%   perround    blocks dispatched per round; 0 = one per worker. Lower it if a
%               round does not fit in RAM: a round holds
%               2 * perround * blocksize * nsamples * 8 bytes on the client,
%               plus about three times a block on each worker.           (0)
%   verbose     report each block of channels as it finishes             (true)
%   sigtest     only subtract a line where it is
%               significant (see below)                  (false)
%   scanbw      width (Hz) of the band searched around each linefreq for the
%               strongest significant peak, instead of testing the nominal bin
%               only. Use when the line wanders (check with a spectrogram
%               around 50 Hz first). 0 disables. Requires sigtest.       (0)
%   iterations  repeat passes, dropping lines per channel once they stop
%               being found or stop improving. Only useful with scanbw:
%               without it, extra passes are measurably a no-op.         (1)
%   pad         FFT padding factor, nfft = 2^(nextpow2(Nwin)+pad).
%               Finer bin spacing lets the fit land closer to the true line
%               frequency: at Fs 500 / 4 s windows, pad 0 gives 0.244 Hz bins
%               and pad 2 gives 0.061 Hz, at ~4x the runtime. Defaults to 0,
%               which is what the legacy path actually runs at -- cleanline.m
%               writes its PaddingFactor to params.g.pad, which getparams()
%               never reads, so the documented default of 2 is silently
%               ignored there.                            (0)
%
% With sigtest = false (the default) the output is numerically identical to
%   pop_cleanline(EEG, 'scanforlines', false, 'computepower', false, ...)
% because with explicit line frequencies fitlinesc() picks its frequency
% bins from the frequency grid alone -- no coupling between channels.
%
% sigtest = true gates the subtraction on the Thomson F-test, per window and
% per channel: a sinusoid at f0 is fitted and removed only where its F-value
% clears the p threshold. The legacy path does not do this -- fitlinesc() takes
% its explicit-f0 branch (fitlinesc.m:112-121) and subtracts at f0
% unconditionally, computing `sig` and never consulting it, so it removes
% energy even from channels and windows that carry no line. This is what the
% newer cleanLineNoise()/PREP implementation does differently
% (fitSignificantFrequencies.m:59).
%
% The gated path also fixes the overlap-smoothing bug at rmlinesmovingwinc.m:106
% (datafitwin0 is captured before the blend rather than after, so each window is
% crossfaded with its own tail instead of the previous window's). That is inert
% while every window fits the same bin, but it breaks scanning outright, and
% without it scanbw made the drift case *worse*. Cost: with sigtest on, output is
% no longer bit-identical to the legacy loop -- about 0.4% of signal RMS on
% static lines, in the direction of slightly deeper removal. The default path
% (sigtest false) still calls rmlinesmovingwinc() untouched and remains exactly
% equal to pop_cleanline().
%
% Measured on 64 ch x 300 s, half the channels carrying 50/100/150 Hz and eight
% carrying a line drifting +-0.4 Hz (50 Hz residual / drifting-line residual /
% RMS disturbance of line-free channels / runtime on 8 workers):
%   no scan                    -24.11  -9.82  0.0127    3 s
%   scanbw 2, iterations 3     -27.44 -24.44  0.0612   15 s
%   scanbw 2, iterations 5     -27.55 -24.90  0.0612   24 s
%   cleanLineNoise reference   -27.55 -24.91  0.0612   32 s
% Note the scan is what disturbs line-free channels (0.0127 -> 0.06): searching
% a 2 Hz band across many windows finds far more spurious significant peaks than
% testing one bin does. Only turn it on if your line actually moves.
%
% See Also: cleanline(), rmlinesmovingwinc(), cleanLineNoise()

% The only external dependency is chronux, in the cleanline plugin's external/
% folder. The caller is expected to have run
%   addpath(genpath('<path to cleanline plugin>'))
% before the parpool is created, so the workers inherit it too.
if ~exist('rmlinesmovingwinc', 'file')
    error('cleanline_fast:noChronux', ...
        ['rmlinesmovingwinc() not found. Add the cleanline plugin to the path ' ...
         'with addpath(genpath(''<path to cleanline plugin>'')).']);
end

g = struct('chanlist',  1:EEG.nbchan, ...
           'linefreqs', 60, ...
           'p',         0.01, ...
           'bandwidth', 2, ...
           'winsize',   4, ...
           'winstep',   1, ...
           'tau',       100, ...
           'blocksize', 4, ...
           'sigtest',   false, ...
           'pad',        0, ...
           'scanbw',     0, ...
           'iterations', 1, ...
           'perround',   0, ...
           'verbose',    true);
for k = 1:2:numel(varargin)
    if ~isfield(g, lower(varargin{k}))
        error('cleanline_fast:badOption', 'Unknown option ''%s''', varargin{k});
    end
    g.(lower(varargin{k})) = varargin{k+1};
end

% Scanning picks the strongest *significant* peak in a band, and iterating an
% ungated subtraction just removes the same estimate repeatedly, so both only
% make sense with the F-test in play.
if ~g.sigtest && (g.scanbw > 0 || g.iterations > 1)
    error('cleanline_fast:needSigtest', ...
        '''scanbw'' and ''iterations'' require ''sigtest'', true.');
end

% A line at or above Nyquist has no bin to fit; min(abs(f-f0)) would silently
% snap it to the top of the grid and subtract at the wrong frequency.
tooHigh = g.linefreqs >= EEG.srate/2;
if any(tooHigh)
    warning('cleanline_fast:aboveNyquist', ...
        'Ignoring line frequencies at or above Nyquist (%g Hz): %s', ...
        EEG.srate/2, mat2str(g.linefreqs(tooHigh)));
    g.linefreqs(tooHigh) = [];
end
if isempty(g.linefreqs)
    error('cleanline_fast:noLines', ...
        'No line frequencies left below Nyquist (%g Hz).', EEG.srate/2);
end

params    = struct('tapers', [g.bandwidth/2, g.winsize, 1], ...
                   'Fs',     EEG.srate, ...
                   'pad',    g.pad);
movingwin = [g.winsize g.winstep];

dataclass = class(EEG.data);

chans  = g.chanlist(:)';
starts = 1:g.blocksize:numel(chans);
blocks = arrayfun(@(i) chans(i:min(i + g.blocksize - 1, numel(chans))), ...
                  starts, 'UniformOutput', false);

% Work in rounds, taking each block straight out of EEG.data and writing it
% straight back, rather than staging the whole recording as one double array.
% For 256 channels x 8 h at 250 Hz that array alone would be ~15 GB, and the
% sliced copies parfor needs would add as much again; this way only one round
% is resident. Trials are concatenated by the (chan, :) indexing, which is how
% the legacy loop treats them too.
% Size the rounds to the pool -- and make sure the pool exists first. gcp with
% 'nocreate' returns empty when none is running yet, which would set perRound to
% 1 and feed the workers one block at a time while the rest idle; parfor would
% then start the pool anyway, so the run looks parallel and is not.
pool = gcp('nocreate');
% A thread pool cannot run this: the multitaper inner loop calls MEX, and thread workers
% reject it ("Use of MEX functions is not supported on a thread-based worker"). The pool
% is session-global, so one left running by another stage - GEDAI's band loop with
% PoolType 'Threads', say - would land here and fail mid-recording, after Zapline has
% already spent its ten minutes. Swap it for a process pool instead.
if ~isempty(pool) && contains(class(pool), 'ThreadPool')
    fprintf('cleanline_fast: thread pool is running; replacing it with a process pool (MEX).\n');
    delete(pool);
    pool = [];
end
if isempty(pool)
    try
        pool = parpool;
        % The next recording spends ~18 min in Zapline with no parfor in sight, which is
        % long enough for the default IdleTimeout to reap this pool; the stage after it
        % then pays the startup again. Repeated 30-worker startup and teardown is also
        % what precedes the failures where only some workers ever connect, so keep the
        % pool alive across recordings instead of rebuilding it for each one.
        try
            pool.IdleTimeout = max(pool.IdleTimeout, 240);
        catch
        end
    catch ME
        % No Parallel Computing Toolbox, or the pool would not come up - asking for the
        % machine's full worker count, parpool can sit at "Connected to 27 of 30 workers"
        % for twenty minutes and then give up. Half a pool beats the serial fallback,
        % which for 250 channels is the difference between minutes and hours.
        fprintf('cleanline_fast: parpool failed (%s); retrying with fewer workers.\n', ...
            ME.message);
        try
            cl   = parcluster;
            pool = parpool(cl, max(4, floor(cl.NumWorkers / 2)));
        catch ME2
            fprintf('cleanline_fast: still no pool (%s); running serially.\n', ME2.message);
            pool = [];
        end
    end
end
if isempty(pool)
    perRound = 1;
else
    perRound = pool.NumWorkers;
end

% parfor's worker limit, NOT the number of workers to ask for. Zero forces the loop to
% run in the client. Without it a bare parfor creates a pool of its own whenever none is
% running, so the catch above buys nothing: the same pool that just failed to start gets
% requested again, and the second attempt queues indefinitely rather than erroring.
if isempty(pool), nWork = 0; else, nWork = pool.NumWorkers; end
if g.perround > 0
    perRound = g.perround;      % lower this if a round does not fit in RAM
end

nb = numel(blocks);
if g.verbose
    fprintf('cleanline_fast: %d channels, %d blocks of %d, %d per round\n', ...
        numel(chans), nb, g.blocksize, perRound);
end
tStart = tic;

for r = 1:perRound:nb
    rb = r:min(r + perRound - 1, nb);

    % One strided gather per round rather than one per block: EEG.data is
    % channels-by-samples, so pulling a few rows walks the whole array.
    chunk = double(EEG.data([blocks{rb}], :)).';
    in    = cell(1, numel(rb));
    c0    = 0;
    for j = 1:numel(rb)
        nc    = numel(blocks{rb(j)});
        in{j} = chunk(:, c0 + 1:c0 + nc);
        c0    = c0 + nc;
    end
    clear chunk;
    out = cell(1, numel(rb));

    bl   = blocks(rb);          % sliced, so the workers can name their channels
    verb = g.verbose;

    parfor (j = 1:numel(rb), nWork)
        x = in{j};
        if g.sigtest
            datac = removelines_gated(x, movingwin, g.tau, params, g.p, ...
                                      g.linefreqs, g.scanbw, g.iterations);
        else
            datac = rmlinesmovingwinc(x, movingwin, g.tau, params, g.p, 'n', g.linefreqs, []);
        end
        % The sliding windows need not tile the recording exactly, so the tail
        % is carried over raw. cleanline.m:405 writes that tail with
        %   datac(end:end+ndiff) = data(end-ndiff:end)
        % whose leading index is off by one, so it also reverts the last
        % *cleaned* sample to its raw value. Reproduced so output matches.
        L = size(datac, 1);
        if L < size(x, 1)
            datac(L:size(x, 1), :) = x(L:end, :);
        end
        out{j} = datac;
        if verb
            % Channels in a block are fitted together, so they finish together.
            % Worker output is forwarded to the client, but the order across
            % workers is whatever finishes first.
            fprintf('Cleaned Chan %s\n', strtrim(sprintf('%d ', bl{j})));
        end
    end

    for j = 1:numel(rb)
        EEG.data(blocks{rb(j)}, :) = cast(out{j}.', dataclass);
    end
    clear in out;

    if g.verbose
        el = toc(tStart);
        fprintf('  %3d/%d blocks | %5.1f min elapsed | ~%4.1f min left\n', ...
            rb(end), nb, el/60, el/60*(nb - rb(end))/rb(end));
    end
end

if g.verbose
    tTotal = toc(tStart);
    if tTotal < 90
        fprintf('cleanline_fast: %d channels done in %.1f s\n', numel(chans), tTotal);
    else
        fprintf('cleanline_fast: %d channels done in %.1f min (%.0f s)\n', ...
            numel(chans), tTotal/60, tTotal);
    end
end

if ~isempty(EEG.icaweights)
    EEG.icaact = [];   % stale now that the channel data changed
end


function datac = removelines_gated(data, movingwin, tau, params, p, f0, scanbw, iters)
% rmlinesmovingwinc() with the subtraction gated on the Thomson F-test,
% optionally scanning a band around each line and repeating the pass.
%
% With scanbw = 0 and iters = 1 this is rmlinesmovingwinc() with exactly one
% thing changed: a sinusoid is fitted only where its F-value clears the p
% threshold for that window and that channel, instead of unconditionally.
%
% scanbw > 0 searches +-scanbw/2 around each line for the strongest peak that
% clears the threshold, which is what catches a line that does not sit still.
% iters > 1 repeats the whole pass; a line is retired for a channel once it is
% no longer found anywhere, or once removing it stopped reducing power at the
% nominal frequency (the over-subtraction guard). Both mirror cleanLineNoise().

data     = change_row_to_column(data);
[N, C]   = size(data);
Fs       = params.Fs;
Nwin     = round(Fs*movingwin(1));
Nstep    = round(movingwin(2)*Fs);
Noverlap = Nwin - Nstep;

xs     = (1:Noverlap)';
smooth = repmat(1./(1 + exp(-tau.*(xs - Noverlap/2)/Noverlap)), [1 C]);

winstart = 1:Nstep:(N - Nwin + 1);
nw       = length(winstart);
Lfit     = winstart(nw) + Nwin - 1;

% tapers and frequency grid are fixed across windows and passes, so build once
[tapers, pad, ~, fpass] = getparams(params);   % resolves [W T p] -> [TW K]
wparams        = params;
wparams.tapers = dpsschk(tapers, Nwin, Fs);
nfft           = max(2^(nextpow2(Nwin) + pad), Nwin);
[f, findx]     = getfgrid(Fs, nfft, fpass);
tt             = (0:Nwin - 1)';

nf   = numel(f0);
fbin = zeros(1, nf);
lo   = zeros(1, nf);
hi   = zeros(1, nf);
for k = 1:nf
    [~, fbin(k)] = min(abs(f - f0(k)));
    if scanbw > 0
        [~, lo(k)] = min(abs(f - (f0(k) - scanbw/2)));
        [~, hi(k)] = min(abs(f - (f0(k) + scanbw/2)));
    end
end

% Bins the convergence guard watches: the whole scanned band when scanning,
% since a line caught at 50.3 Hz need not reduce power at the 50.0 Hz bin.
probe = cell(1, nf);
for k = 1:nf
    if scanbw > 0
        probe{k} = lo(k):hi(k);
    else
        probe{k} = fbin(k);
    end
end

datac  = data;
active = true(nf, C);        % lines still worth chasing, per channel
if iters > 1
    Pprev = bandpower_db(datac, Nwin, wparams, Fs, nfft, findx, probe);
end

for it = 1:iters
    datafit = zeros(Lfit, C);
    found   = false(nf, C);  % line was significant somewhere during this pass

    for n = 1:nw
        idx = winstart(n):(winstart(n) + Nwin - 1);
        [Fval, A, ~, sig] = ftestc(datac(idx, :), wparams, p, 'n');

        datafitwin = zeros(Nwin, C);
        for ch = 1:C
            ka = find(active(:, ch)).';
            if isempty(ka)
                continue;
            end
            bins = zeros(1, numel(ka));
            nb   = 0;
            for k = ka
                if scanbw > 0
                    r     = lo(k):hi(k);
                    Fscan = Fval(r, ch);
                    Fscan(Fscan < sig) = 0;
                    if ~any(Fscan)
                        continue;            % nothing significant in the band
                    end
                    [~, im] = max(Fscan);
                    b = r(1) + im - 1;
                else
                    if Fval(fbin(k), ch) < sig
                        continue;            % <-- the F-test gate
                    end
                    b = fbin(k);
                end
                nb          = nb + 1;
                bins(nb)    = b;
                found(k, ch) = true;
            end
            if nb == 0
                continue;                    % no line here: subtract nothing
            end
            bins = unique(bins(1:nb), 'stable');   % overlapping bands may collide
            fsig = reshape(f(bins), 1, []);
            datafitwin(:, ch) = real(exp(1i*2*pi*tt*fsig/Fs)*A(bins, ch) + ...
                                     exp(-1i*2*pi*tt*fsig/Fs)*conj(A(bins, ch)));
        end

        % Overlap smoothing. rmlinesmovingwinc.m:106 captures datafitwin0 before
        % the blend rather than after, so it crossfades each window with its own
        % tail instead of the previous window's -- inert when every window fits
        % the same bin, but actively wrong once scanning lets adjacent windows
        % fit different frequencies. Done properly here.
        prev = datafitwin;
        if n > 1
            datafitwin(1:Noverlap, :) = smooth.*datafitwin(1:Noverlap, :) + ...
                (1 - smooth).*prev0((Nwin - Noverlap + 1):Nwin, :);
        end
        prev0 = prev; %#ok<NASGU>
        datafit(idx, :) = datafitwin;
    end

    datac(1:Lfit, :) = datac(1:Lfit, :) - datafit;

    if iters > 1
        Pnow   = bandpower_db(datac, Nwin, wparams, Fs, nfft, findx, probe);
        worse  = (Pprev - Pnow) < 0;      % this pass did not reduce the line
        Pprev  = Pnow;
        active = active & found & ~worse;
        if ~any(active(:))
            break;
        end
    end
end

datac = datac(1:Lfit, :);


function P = bandpower_db(x, Nseg, wparams, Fs, nfft, findx, probe)
% Peak multitaper power (dB) within each probe band, averaged over
% non-overlapping Nseg-sample segments. Batched across channels, unlike
% calculateSegmentSpectrum(), which is univariate.
[N, C] = size(x);
nf     = numel(probe);
ns     = floor(N/Nseg);
if ns < 1
    P = -inf(nf, C);
    return;
end
X = reshape(x(1:ns*Nseg, :), Nseg, ns*C);
J = mtfftc(X, wparams.tapers, nfft, Fs);
J = J(findx, :, :);
S = mean(conj(J).*J, 2);                  % average over tapers
S = reshape(S, numel(findx), ns, C);
S = reshape(mean(S, 2), numel(findx), C); % average over segments
S = abs(S);
P = zeros(nf, C);
for k = 1:nf
    P(k, :) = 10*log10(max(S(probe{k}, :), [], 1));
end
