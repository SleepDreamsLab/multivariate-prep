function [failures, GEDs] = bidsfun_subcomp(BIDS, opts)
% BIDSFUN_SUBCOMP  Subtract the ICLabel-flagged ICA components, then run a GED.
%
%   Per recording: load the bidsfun_gedai .set and the run-pamica.py _ica.mat,
%   attach the decomposition (loadica), subtract the components flagged as artefact
%   with pop_subcomp, and hand the cleaned recording to ged() - the generalized
%   eigendecomposition (see GED/ged.m). Nothing is written to disk: the GED results
%   are returned so the caller can inspect them and decide what is worth keeping.
%
%   Reads   <inputpath>/<sub>/<ses>/<fileID>_desc-<inputdesc>_eeg<inputfileext>
%           <icapath>/<sub>/<ses>/<fileID>_desc-<icadesc>_ica.mat   (+ _iclabels.tsv)
%   Writes  nothing.
%
%   The flags come from the _iclabels.tsv `status` column when preferstatus is true
%   (the default - that is the file meant to be hand-screened), otherwise from the
%   .mat's gcompreject; see loadica. A missing input is "not ready yet" (skipped,
%   not an error), so run-pamica.py can still be writing it on another machine.
%
%   Name-value, all optional:
%     derivfolder                     derivatives subfolder for the inputs ('prep-ged')
%     inputpath/inputdesc/inputfileext   the EEG          ('' -> deriv, 'filt2ged', '.set')
%     icapath/icadesc                     the _ica.mat    ('' -> inputpath, 'pamica')
%     preferstatus                    obey the .tsv status column over gcompreject   (true)
%     gedargs                         name-value cell forwarded to ged(), e.g.
%                                     {'contrast','spectral','peakfreq',13.5,'fwhm',2}
%                                     ({} -> ged's own defaults: a 10 Hz spectral
%                                     contrast against broadband)
%     tasklabel/acqlabel/subjectfilter/sessionfilter        BIDS query and filters
%
%   Outputs:
%     failures   cell of structs, one per recording that errored
%     GEDs       struct array with fields fileID, removed (the subtracted ICs) and
%                ged (the struct ged() returned)
%
% Methods section:
%
% The independent components classified as artefact (see bidsfun_iclabel /
% run-pamica.py) were removed by subtracting their back-projection from the sensor
% data. The cleaned recordings then entered a generalized eigendecomposition
% (Cohen, 2022) to derive spatial filters maximising the contrast of interest.

arguments
    BIDS
    opts.derivfolder   char = 'prep-ged'
    opts.inputpath     char = ''
    opts.inputdesc     char = 'filt2ged'
    opts.inputfileext  char = '.set'
    opts.icapath       char = ''
    opts.icadesc       char = 'pamica'
    opts.preferstatus (1,1) logical = true
    opts.gedargs       cell = {}
    opts.tasklabel          = {'Sleep', 'sleep'}
    opts.acqlabel      char = ''
    opts.subjectfilter cell = {}
    opts.sessionfilter cell = {}
end

fprintf('\n=== Running bidsfun_subcomp ===\n');

if isempty(opts.inputpath), opts.inputpath = fullfile(BIDS.pth, 'derivatives', opts.derivfolder); end
if isempty(opts.icapath),   opts.icapath   = opts.inputpath; end

filesEEG = bids.query(BIDS, 'data', 'extension', '.vhdr', ...
    'task', opts.tasklabel, 'acq', opts.acqlabel);
if isempty(filesEEG)
    error('bidsfun_subcomp:noFiles', 'No matching EEG files found in BIDS layout.');
end

failures = {};
GEDs     = struct('fileID', {}, 'removed', {}, 'ged', {});
for ifile = 1:numel(filesEEG)
    p      = bids.internal.parse_filename(filesEEG{ifile});
    fileID = strjoin(cellfun(@(k) [k '-' p.entities.(k)], fieldnames(p.entities), 'uni', 0), '_');
    subDir = fullfile(['sub-' p.entities.sub], ['ses-' p.entities.ses]);

    if ~isempty(opts.subjectfilter) && ~contains(fileID, opts.subjectfilter), continue, end
    if ~isempty(opts.sessionfilter) && ~contains(fileID, opts.sessionfilter), continue, end
    fprintf('\n=== %s ===\n', fileID)

    inFile  = fullfile(opts.inputpath, subDir, [fileID '_desc-' opts.inputdesc '_eeg' opts.inputfileext]);
    icaFile = fullfile(opts.icapath,   subDir, [fileID '_desc-' opts.icadesc   '_ica.mat']);

    if ~isfile(inFile) || ~isfile(icaFile)
        fprintf('[skip] input not ready (EEG %d, ICA %d)\n', isfile(inFile), isfile(icaFile)); continue
    end

    try
        EEG = fast_eeg_import(inFile);
        EEG = loadica(EEG, icaFile, 'preferstatus', opts.preferstatus, 'checkset', false);

        if ~isfield(EEG.reject, 'gcompreject') || isempty(EEG.reject.gcompreject)
            error('bidsfun_subcomp:noFlags', ...
                '%s has no artefact flags - run ICLabel first.', icaFile);
        end
        badComps = find(EEG.reject.gcompreject);
        fprintf('Subtracting %d/%d components: %s\n', ...
            numel(badComps), size(EEG.icaweights, 1), mat2str(badComps));
        if ~isempty(badComps)
            EEG = pop_subcomp(EEG, badComps, 0);
        end

        %%% The decomposition no longer describes the data - drop it so nothing
        %%% downstream subtracts a second time; keep a record of what was removed.
        EEG.etc.ic_subtraction = struct('icaFile', icaFile, 'removed', badComps);
        [EEG.icaweights, EEG.icasphere, EEG.icawinv, EEG.icachansind, EEG.icaact] = deal([]);

        %%% GED on the cleaned recording. Which contrast to run is the caller's
        %%% decision (gedargs) - it is the choice that determines what the spatial
        %%% filters end up isolating.
        GED = ged(EEG, opts.gedargs{:});
        GEDs(end+1) = struct('fileID', fileID, 'removed', badComps, 'ged', GED); %#ok<AGROW>

    catch ME
        fprintf('[ERROR] %s: %s\n', fileID, ME.message);
        failures{end+1} = struct('fileID', fileID, 'message', ME.message); %#ok<AGROW>
    end
end

if ~isempty(failures)
    fprintf('\n=== %d file(s) failed ===\n', numel(failures));
    for k = 1:numel(failures)
        fprintf('  %s: %s\n', failures{k}.fileID, failures{k}.message);
    end
end
end
