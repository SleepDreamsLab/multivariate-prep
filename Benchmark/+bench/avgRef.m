function data = avgRef(data)
% AVGREF  Average reference, counting the recording reference as an implicit zero channel.
%
%   data = bench.avgRef(data)
%
%   The same operation bidsfun_gedai applies before GEDAI: the reference electrode (Cz on
%   EGI nets) is not a row of the data but is part of the montage, so the mean runs over
%   nChan + 1 channels. The operation is linear, so avgRef(B + S) == avgRef(B) + avgRef(S):
%   background and injected ground truth can be referenced separately and still add up
%   to exactly what GEDAI sees.
%
%   data   [nChan x nSamples], referenced to the recording reference. Returned re-referenced.

data = data - sum(data, 1) / (size(data, 1) + 1);
end
