function [freeBytes, totalBytes] = freemem()
% FREEMEM  Available and total physical RAM, in bytes, on Windows/Linux/macOS.
%
%   [freeBytes, totalBytes] = utils.freemem()
%
%   MATLAB's built-in MEMORY is Windows-only ("Function MEMORY is not available on
%   this platform." on Linux/macOS), so anything that sizes a workload against free
%   RAM has to ask the OS itself off Windows.
%
%   Linux : /proc/meminfo MemAvailable - the kernel's own estimate of what a new
%           allocation can get without swapping (MemFree plus the reclaimable part of
%           page cache and slab). MemFree alone would badly understate it, since Linux
%           parks nearly all idle RAM in page cache.
%   macOS : vm_stat - free + inactive + speculative + purgeable pages, the closest
%           analogue to MemAvailable; total from sysctl hw.memsize.
%
%   Either value is NaN if the platform is unrecognised or the query fails; callers
%   must handle that rather than silently sizing off a bogus number.
%
%   Note: on Linux this reports the *host's* memory. Inside a container or under a
%   cgroup/SLURM memory limit the real ceiling can be lower, and this will overstate
%   what is available.

freeBytes  = NaN;
totalBytes = NaN;

try
    if ispc
        [~, sysMem] = memory;
        freeBytes   = sysMem.PhysicalMemory.Available;
        totalBytes  = sysMem.PhysicalMemory.Total;

    elseif ismac
        [st, out] = system('sysctl -n hw.memsize');
        if st == 0, totalBytes = str2double(strtrim(out)); end

        [st, out] = system('vm_stat');
        if st == 0
            pageSize = regexp(out, 'page size of (\d+) bytes', 'tokens', 'once');
            pageSize = str2double(pageSize{1});
            nPages   = @(name) local_pages(out, name);
            freeBytes = pageSize * (nPages('free') + nPages('inactive') + ...
                nPages('speculative') + nPages('purgeable'));
        end

    elseif isunix
        txt = fileread('/proc/meminfo');
        freeBytes  = local_meminfo(txt, 'MemAvailable');
        totalBytes = local_meminfo(txt, 'MemTotal');
        if isnan(freeBytes)
            %%% Kernels before 3.14 have no MemAvailable field; MemFree + Cached is the
            %%% usual stand-in (coarser, since not all page cache is reclaimable).
            freeBytes = local_meminfo(txt, 'MemFree') + local_meminfo(txt, 'Cached');
        end
    end
catch
    % leave NaN - the caller decides what to do without a reading
end
end

% -------------------------------------------------------------------------
function bytes = local_meminfo(txt, key)
% /proc/meminfo lines look like "MemAvailable:   12345678 kB".
tok = regexp(txt, ['^' key ':\s*(\d+)\s*kB'], 'tokens', 'once', 'lineanchors');
if isempty(tok), bytes = NaN; else, bytes = str2double(tok{1}) * 1024; end
end

% -------------------------------------------------------------------------
function n = local_pages(txt, name)
% vm_stat lines look like "Pages free:                          123456."
tok = regexp(txt, ['Pages ' name ':\s*(\d+)\.'], 'tokens', 'once');
if isempty(tok), n = 0; else, n = str2double(tok{1}); end
end
