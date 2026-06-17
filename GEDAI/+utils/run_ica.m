function EEG = run_ica(EEG, opts)
% GEDAI.UTILS.RUN_ICA  Rank-aware ICA on a continuous EEGLAB structure.
%
%   EEG = gedai.utils.run_ica(EEG)
%   EEG = gedai.utils.run_ica(EEG, Name, Value)
%
%   Runs AMICA (preferred) or extended infomax on EEG.data.  The effective
%   rank of the data is determined automatically from the covariance
%   eigenvalue spectrum, so the function handles both full-rank data and
%   rank-deficient data (e.g. after GEDAI projects out GED components).
%   Call this function twice — once for the raw-stage data, once for the
%   GEDAI-cleaned data — with whatever EEG structure you want decomposed.
%
%   Inputs
%   ------
%   EEG              : Continuous EEGLAB EEG structure (EEG.trials == 1).
%
%   Optional name-value pairs
%   -------------------------
%   Rank             : Integer. Effective rank to use for ICA.
%                      Default: auto-detected from the eigenvalue cliff.
%   RankThreshold    : Relative eigenvalue threshold for auto-detection
%                      (default: 1e-6, i.e. eigenvalues < threshold *
%                      max_eigenvalue are considered zero).
%   Method           : 'amica' | 'infomax' (default: 'amica' if runamica
%                      is on the path, otherwise 'infomax').
%   AmicaOutDir      : Directory for AMICA output files
%                      (default: tempname in the system temp folder).
%   AmicaNumProcs    : Number of parallel processes for AMICA (default: 4).
%   Verbose          : true | false — print progress (default: true).
%
%   Output
%   ------
%   EEG : Input structure with ICA fields populated:
%         EEG.icaweights, EEG.icasphere, EEG.icawinv, EEG.icaact.
%
%   Sanity check
%   ------------
%   Before calling, you can verify the rank yourself:
%       cov_eig = eig(cov(EEG.data'));
%       plot(sort(cov_eig, 'descend'));  % cliff shows effective rank
%       % or simply:
%       disp(rank(EEG.data))
%
%   Example
%   -------
%   [EEGclean, EEGstage, epochIdx] = gedai.stage_specific(EEG, Scoring, [0]);
%   EEGclean = gedai.utils.run_ica(EEGclean);          % after GEDAI
%   EEGstage = gedai.utils.run_ica(EEGstage);          % before GEDAI

arguments
    EEG
    opts.Rank             (1,1) double = NaN     % NaN = auto-detect
    opts.RankThreshold    (1,1) double = 1e-6
    opts.Method           (1,1) string = "auto"  % 'amica'|'infomax'|'auto'
    opts.AmicaOutDir      (1,:) char   = ''
    opts.AmicaNumProcs    (1,1) double = 4
    opts.Verbose          (1,1) logical = true
end

%%% ------------------------------------------------------------------ %%%
%%% 1.  Detect effective rank                                           %%%
%%% ------------------------------------------------------------------ %%%
if isnan(opts.Rank)
    covMatrix  = cov(double(EEG.data)');
    eigVals    = sort(eig(covMatrix), 'descend');
    eigVals    = eigVals / eigVals(1);                    % normalise
    effectiveRank = sum(eigVals > opts.RankThreshold);

    if opts.Verbose
        fprintf('gedai.utils.run_ica: auto-detected rank = %d  (out of %d channels)\n', ...
            effectiveRank, EEG.nbchan);
        fprintf('  Tip: plot(sort(eig(cov(EEG.data'')),''descend'')) to inspect the cliff.\n');
    end
else
    effectiveRank = opts.Rank;
    if opts.Verbose
        fprintf('gedai.utils.run_ica: using user-supplied rank = %d\n', effectiveRank);
    end
end

if effectiveRank < 2
    error('gedai:utils:run_ica:rankTooLow', ...
        'Effective rank is %d — check your data.', effectiveRank);
end

%%% ------------------------------------------------------------------ %%%
%%% 2.  Choose ICA method                                               %%%
%%% ------------------------------------------------------------------ %%%
if opts.Method == "auto"
    if exist('runamica15', 'file') || exist('runamica', 'file')
        method = "amica";
    else
        method = "infomax";
        if opts.Verbose
            warning('gedai:utils:run_ica:noAmica', ...
                'AMICA (runamica15) not found on path — falling back to extended infomax.');
        end
    end
else
    method = opts.Method;
end

%%% ------------------------------------------------------------------ %%%
%%% 3.  Run ICA                                                         %%%
%%% ------------------------------------------------------------------ %%%
switch method

    % ----------------------------------------------------------------- %
    case "amica"
        outDir = opts.AmicaOutDir;
        if isempty(outDir)
            outDir = fullfile(tempdir, ['amica_' datestr(now, 'yyyymmdd_HHMMSS')]); %#ok<TNOW1,DATST>
        end
        if ~exist(outDir, 'dir'), mkdir(outDir); end

        if opts.Verbose
            fprintf('gedai.utils.run_ica: running AMICA  (rank=%d, procs=%d)\n  outDir: %s\n', ...
                effectiveRank, opts.AmicaNumProcs, outDir);
        end

        amicaFn = 'runamica15';
        if ~exist(amicaFn, 'file'), amicaFn = 'runamica'; end

        feval(amicaFn, double(EEG.data), ...
            'num_chans',  EEG.nbchan, ...
            'num_frames', EEG.pnts, ...
            'outdir',     outDir, ...
            'numprocs',   opts.AmicaNumProcs, ...
            'max_threads', opts.AmicaNumProcs, ...
            'num_mix',    1, ...            % single mixture model
            'do_opt_in',  1, ...
            'pcakeep',    effectiveRank);   % <- rank reduction

        if opts.Verbose
            fprintf('gedai.utils.run_ica: loading AMICA results from %s\n', outDir);
        end
        modout = loadmodout15(outDir);
        EEG.icaweights = modout.W;
        EEG.icasphere  = modout.S;

    % ----------------------------------------------------------------- %
    case "infomax"
        if opts.Verbose
            fprintf('gedai.utils.run_ica: running extended infomax  (rank=%d)\n', effectiveRank);
        end

        [weights, sphere] = runica(double(EEG.data), ...
            'extended', 1, ...
            'pca',      effectiveRank, ...
            'verbose',  char(string(opts.Verbose)));
        EEG.icaweights = weights;
        EEG.icasphere  = sphere;

    otherwise
        error('gedai:utils:run_ica:unknownMethod', ...
            'Unknown Method ''%s''. Use ''amica'' or ''infomax''.', method);
end

%%% ------------------------------------------------------------------ %%%
%%% 4.  Derive inverse matrix and activations                          %%%
%%% ------------------------------------------------------------------ %%%
EEG.icawinv = pinv(EEG.icaweights * EEG.icasphere);
EEG         = eeg_checkset(EEG);

if opts.Verbose
    fprintf('gedai.utils.run_ica: done. %d components on %d channels.\n', ...
        size(EEG.icaweights, 1), EEG.nbchan);
end

end
