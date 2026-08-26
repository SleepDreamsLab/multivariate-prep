function EEG = extendICAToDerived(EEG, EEGpre, chanOpts, derivedIdx)
% EXTENDICATODERIVED  Give chans1020's synthetic EOG channels their own rows
% in the ICA mixing matrix, so that projecting components out cleans them
% too.
%
%   EEG = viewer.extendICAToDerived(EEG, EEGpre, chanOpts, derivedIdx)
%
%   EEG         the recording AFTER chans1020 added the derived channels
%   EEGpre      the same recording BEFORE that call (carries the ICA fields
%               and the original 'E<n>' labels chans1020 looks channels up by)
%   chanOpts    the exact argument list that was passed to chans1020, so the
%               probe below reproduces the same derivation
%   derivedIdx  channel indices of the derived channels (e.g. chanmap.EOG1/2)
%
% WHY THIS EXISTS
% chans1020 appends EOG1/EOG2 as bipolar differences of net channels. Those
% channels are not part of the decomposition's channel set, so a viewer that
% cleans by projecting components out of icachansind leaves them completely
% untouched -- the eye components stay fully visible in precisely the two
% traces you look at to judge whether the eye components were removed.
%
% HOW
% A derived channel is a LINEAR combination of net channels (a bipolar
% difference, whose anchors may themselves be spherical-spline interpolated
% from surviving channels -- also linear). The component map of a linear
% combination of channels is the same linear combination of those channels'
% component maps. So instead of reimplementing chans1020's derivation here
% -- and having to keep it in sync forever, including the interpolation
% branch -- push the MIXING MATRIX through the very same chans1020 call,
% treating one component map as one "sample". The rows that come back for
% the derived channels are exactly their mixing rows.
%
% The unmixing gets zero columns for the derived channels: they are
% redundant with the net channels they are built from, so they must not
% contribute to the activations (that would double-count them). Only the
% mixing side grows, which is precisely what cleaning needs.
if ~isfield(EEG, 'icawinv') || isempty(EEG.icawinv) || isempty(derivedIdx)
    return
end
derivedIdx = derivedIdx(:).';
nComp      = size(EEG.icawinv, 2);
nDer       = numel(derivedIdx);

% Probe dataset: component maps as data, on the pre-chans1020 channel set.
% Channels outside the decomposition get a zero map, which matches how they
% are treated everywhere else -- they simply are not cleaned.
probe = EEGpre;
for f = intersect(fieldnames(probe), {'icaweights','icasphere','icawinv','icachansind','icaact'})'
    probe = rmfield(probe, f{1});
end
probe.data = zeros(EEGpre.nbchan, nComp);
probe.data(EEGpre.icachansind, :) = EEG.icawinv;
probe.pnts   = nComp;
probe.trials = 1;
probe.xmin   = 0;
probe.xmax   = max(nComp - 1, 0) / probe.srate;
probe.event  = [];
probe.epoch  = [];

probeOut = chans1020(probe, chanOpts{:});
if size(probeOut.data, 1) < max(derivedIdx)
    warning('extendICAToDerived:noDerivedRows', ...
        'The chans1020 probe returned %d channels, fewer than the requested derived index %d; leaving the decomposition unextended.', ...
        size(probeOut.data, 1), max(derivedIdx));
    return
end
winvDerived = double(probeOut.data(derivedIdx, :));    % nDer x nComp

EEG.icawinv     = [EEG.icawinv;   winvDerived];
EEG.icasphere   = [EEG.icasphere, zeros(size(EEG.icasphere, 1), nDer)];
EEG.icachansind = [EEG.icachansind, derivedIdx];
% NOTE: deliberately no eeg_checkset(EEG,'ica') after this point --
% checkset recomputes icawinv as pinv(icaweights*icasphere), and the zero
% unmixing columns would send the rows we just added straight back to zero.
end
