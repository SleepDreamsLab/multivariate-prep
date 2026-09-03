function idx = reflectidx(idx, n)
% REFLECTIDX  Fold sample indices back into 1:n by mirroring at both ends.
%
%   Shared by tfmorlet and tfmulti, which both have to read a window of context
%   around a sample and would otherwise run off the ends of the recording.
%   Mirroring rather than zero-padding matters: a taper or a wavelet convolved
%   against a step down to zero produces broadband edge artefacts that look like
%   real events.
%
%   Indexing modulo 2n and folding the upper half back down does it in one
%   vectorised pass over whatever shape it is handed, so the padding costs
%   nothing next to the FFT it feeds. Index 0 maps to 1, -1 to 2, n+1 to n, and
%   so on - whole-sample symmetry, the same convention as 'symmetric' padding.

idx  = mod(idx - 1, 2 * n);
fold = idx >= n;
idx(fold) = 2 * n - 1 - idx(fold);
idx = idx + 1;
end
