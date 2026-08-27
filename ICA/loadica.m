function EEG = loadica(EEG, icafile, opts)
% LOADICA  Attach a run-pamica.py / bidsfun_iclabel decomposition to an EEG struct.
%
%   EEG = loadica(EEG, icafile)
%   EEG = loadica(EEG, icafile, Name, Value, ...)
%
%   Takes the <fileID>_desc-<desc>_ica.mat written by run-pamica.py (optionally
%   already carrying ICLabel results, whether they came from run-pamica.py's
%   mne-icalabel pass or from bidsfun_iclabel.m) and populates the EEGLAB fields
%   on EEG, so the next line can be pop_subcomp / pop_selectcomps / pop_viewprops
%   and everything downstream behaves like a native EEGLAB decomposition.
%
%   Sets, always:
%     EEG.icaweights, EEG.icasphere, EEG.icawinv, EEG.icachansind
%   Sets, when the .mat carries them:
%     EEG.etc.ic_classification.ICLabel   .classes / .classifications / .version
%     EEG.reject.gcompreject              1 x ncomp logical
%     EEG.etc.amica.varfrac               share of back-projected variance per IC
%
%   Typical use:
%     EEG = fast_eeg_import(setFile);
%     EEG = loadica(EEG, icaFile);
%     EEG = pop_subcomp(EEG, find(EEG.reject.gcompreject), 0);
%
%   Required
%   --------
%   EEG       EEGLAB struct, already loaded, with chanlocs.
%   icafile   Path to the _ica.mat. The matching _iclabels.tsv, if present beside
%             it, is used as the fallback source of the artefact flags (see
%             preferstatus below).
%
%   Name-Value
%   ----------
%   preferstatus  Take the artefact flags from the .tsv's `status` column rather
%                 than from the .mat's gcompreject. The .tsv is the file meant to
%                 be hand-edited when a call needs overruling, so this is the
%                 option to use once someone has screened the components.
%                 Default: false (use gcompreject, fall back to the .tsv when the
%                 .mat has no flags).
%   checkset      Run eeg_checkset(EEG, 'ica') at the end. Default: false --
%                 that computes EEG.icaact for the whole recording, ~7 GB for a
%                 236-component night, and nothing needs it: pop_subcomp and
%                 eeg_getdatact derive the activations they want on the fly.
%                 Turn it on only if you actually want icaact resident.
%   verbose       Print what was attached. Default: true.

arguments
    EEG      struct
    icafile  {mustBeTextScalar}
    opts.preferstatus (1,1) logical = false
    opts.checkset     (1,1) logical = false
    opts.verbose      (1,1) logical = true
end

icafile = char(icafile);
if ~isfile(icafile)
    error('loadica:noFile', 'No such file: %s', icafile);
end
S = load(icafile);

for f = {'icaweights', 'icasphere', 'icawinv', 'chanlabels'}
    if ~isfield(S, f{1})
        error('loadica:badFile', ...
            '%s has no "%s" field - is it a run-pamica.py _ica.mat?', icafile, f{1});
    end
end

icalabs = string(S.chanlabels(:));
nComp   = size(S.icaweights, 1);
if size(S.icasphere, 2) ~= numel(icalabs)
    error('loadica:sphereLabelMismatch', ...
        'icasphere has %d columns but chanlabels lists %d channels.', ...
        size(S.icasphere, 2), numel(icalabs));
end

%%% Channel matching. Strict: the unmixing icaweights*icasphere needs every column
%%% it was fitted with, and re-estimating the activations from a subset of channels
%%% by least squares gives something that is not this decomposition.
if isempty(EEG.chanlocs)
    error('loadica:noChanlocs', 'EEG has no chanlocs, so the ICA channels cannot be matched.');
end
eeglabs   = string({EEG.chanlocs.labels}');
[tf, idx] = ismember(icalabs, eeglabs);
if ~all(tf)
    error('loadica:missingChannels', ...
        ['%d of %d channels the decomposition needs are absent from this EEG (%s). ' ...
         'Either the ICA and the EEG are not the same recording, or channels were ' ...
         'dropped after AMICA ran.'], ...
        nnz(~tf), numel(tf), strjoin(icalabs(~tf), ', '));
end

if isfield(S, 'srate') && ~isempty(S.srate) && abs(double(S.srate) - EEG.srate) > 1e-6
    warning('loadica:srateMismatch', ...
        'The decomposition was fitted at %g Hz but this EEG is sampled at %g Hz.', ...
        double(S.srate), EEG.srate);
end

EEG.icaweights  = double(S.icaweights);
EEG.icasphere   = double(S.icasphere);
EEG.icawinv     = double(S.icawinv);
EEG.icachansind = idx(:)';
EEG.icaact      = [];

%%% ICLabel, in the field layout EEGLAB's own iclabel() uses - which is what makes
%%% pop_selectcomps and pop_viewprops title each component with its class.
hasLabels = isfield(S, 'ic_classification');
if hasLabels
    EEG.etc.ic_classification = S.ic_classification;
    cls = S.ic_classification.ICLabel;
    if size(cls.classifications, 1) ~= nComp
        error('loadica:classificationSize', ...
            'ic_classification holds %d rows for %d components.', ...
            size(cls.classifications, 1), nComp);
    end
end

%%% Artefact flags. gcompreject is what the labelling stage decided; the .tsv's
%%% status column is what a human may since have overruled it with, so that one
%%% wins when asked for.
tsvfile = regexprep(icafile, '_ica\.mat$', '_iclabels.tsv');
flags   = [];
source  = '';
if opts.preferstatus && isfile(tsvfile)
    flags = readStatus(tsvfile, nComp);
    source = 'iclabels.tsv (status column)';
elseif isfield(S, 'gcompreject') && ~isempty(S.gcompreject)
    flags = logical(S.gcompreject(:))';
    source = '_ica.mat (gcompreject)';
elseif isfile(tsvfile)
    flags = readStatus(tsvfile, nComp);
    source = 'iclabels.tsv (status column)';
end
if ~isempty(flags)
    EEG.reject.gcompreject = flags;
end

if isfield(S, 'varfrac') && ~isempty(S.varfrac)
    EEG.etc.amica.varfrac = double(S.varfrac(:));
end
EEG.etc.amica.icafile = icafile;

if opts.checkset
    EEG = eeg_checkset(EEG, 'ica');
end

if opts.verbose
    fprintf('loadica: %d components over %d channels from %s\n', ...
        nComp, numel(icalabs), icafile);
    if hasLabels
        probs = S.ic_classification.ICLabel.classifications;
        names = S.ic_classification.ICLabel.classes;
        [~, top] = max(probs, [], 2);
        fprintf('  ICLabel (%s):', S.ic_classification.ICLabel.version);
        for c = 1:numel(names)
            if any(top == c), fprintf(' %s=%d', names{c}, nnz(top == c)); end
        end
        fprintf('\n');
    else
        fprintf('  no ICLabel results in this file\n');
    end
    if ~isempty(flags)
        fprintf('  %d/%d components flagged bad, from %s\n', nnz(flags), nComp, source);
    end
end
end

% -------------------------------------------------------------------------
function flags = readStatus(tsvfile, nComp)
% Read the status column out of an _iclabels.tsv as a 1 x nComp logical.
T = readtable(tsvfile, 'FileType', 'text', 'Delimiter', '\t', ...
    'VariableNamingRule', 'preserve');
if ~ismember('status', T.Properties.VariableNames)
    error('loadica:noStatusColumn', '%s has no "status" column.', tsvfile);
end
if height(T) ~= nComp
    error('loadica:tsvSize', ...
        '%s has %d rows but the decomposition has %d components.', ...
        tsvfile, height(T), nComp);
end
status = string(T.status);
bad    = status == "bad";
if ~all(bad | status == "good")
    error('loadica:badStatusValue', ...
        '%s: the status column must be "good" or "bad"; found %s.', ...
        tsvfile, strjoin(unique(status(~(bad | status == "good"))), ', '));
end
flags = reshape(bad, 1, []);
end
