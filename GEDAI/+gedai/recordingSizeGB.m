function gb = recordingSizeGB(vhdrFile)
% RECORDINGSIZEGB  Size of a BrainVision recording once loaded into memory, in GB.
%   gb = gedai.recordingSizeGB(vhdrFile)
%
%   Derived from the header and the binary's size on disk - no sample data is read, so
%   this is cheap enough to call on a network share before a batch starts.
%
%   EEGLAB holds samples as single, so a 16-bit file doubles on load and the number that
%   matters is not the size on disk. When the header does not name a format the narrowest
%   one is assumed, which errs towards a larger estimate - and a larger estimate means a
%   smaller, safer parallel pool.
%
%   Returns NaN when the header or the binary cannot be read. Callers should treat that
%   as 'unknown' and fall back to a memory-only rule, never as zero.
gb = NaN;
try
    txt      = fileread(vhdrFile);
    dataName = regexp(txt, 'DataFile\s*=\s*([^\r\n]+)', 'tokens', 'once');
    fmt      = regexp(txt, 'BinaryFormat\s*=\s*([^\r\n]+)', 'tokens', 'once');
    if isempty(dataName), return; end
    d = dir(fullfile(fileparts(vhdrFile), strtrim(dataName{1})));
    if isempty(d), return; end
    onDisk = 2;
    if ~isempty(fmt)
        switch upper(strtrim(fmt{1}))
            case {'INT_32', 'IEEE_FLOAT_32'}, onDisk = 4;
            otherwise,                        onDisk = 2;
        end
    end
    gb = (d.bytes * (4 / onDisk)) / 2^30;
catch
end
end
