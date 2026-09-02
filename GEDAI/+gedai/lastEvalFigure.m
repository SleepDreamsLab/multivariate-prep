function suffix = lastEvalFigure(plotflags)
% LASTEVALFIGURE  Filename suffix of the last figure run.eval_clean will write.
%   suffix = gedai.lastEvalFigure(plotflags)
%
%   plotflags is any struct carrying the Plot* toggles - bidsfun_evalfigs' opts, or the
%   name-value struct handed to run.eval_clean. Absent fields count as disabled. Returns
%   '' when no figure with a predictable name is requested.
%
%   This exists so a resume check can ask "did this recording finish?" instead of "did it
%   start?". Keying on a fixed name only answers the second question: whichever figure is
%   written first lands long before the run is over, so a recording that died halfway
%   through looks identical to a complete one and is skipped for good.
%
%   The table below is the write order in run.eval_clean - keep the two in step when a
%   plot is added, moved or renamed.
%
%   Two entries are marked unpredictable and are never chosen as the sentinel:
%
%   PlotCharacteristics  writes only when the cleaned recording carries an etc.GEDAI (or
%                        etc.gedai) field, and eval_clean swallows the error when it does
%                        not. For a pre-GEDAI stage the toggle is on and the file never
%                        appears, so waiting on it would mean never skipping anything.
%   PlotEpochOverlay     writes one file per plotted epoch, named for the epoch index and
%                        its stage label - which depend on the scoring, so the name is not
%                        knowable before the recording has been read.
%
%   When one of those is the last enabled plot, the last predictable figure before it is
%   used instead. A run that then died inside that final step still counts as complete;
%   the alternative, never skipping anything, is worse.

arguments
    plotflags struct
end

%%% flag name              suffix written             predictable name?
order = { ...
    'PlotCharacteristics',  'gedai_characteristics',  false ; ...
    'PlotPsdPerStage',      'psd_per_stage',          true  ; ...
    'PlotPsdPerStageChans', 'psd_per_stage_chans',    true  ; ...
    'PlotPsdOverview',      'psd_overview',           true  ; ...
    'PlotTopoBandPower',    'topo_band_power',        true  ; ...
    'PlotTopoBandStage',    'topo_band_stage_clean',  true  ; ...
    'PlotEpochOverlay',     '',                       false ; ...
    'PlotTimefreq',         'fooof2-45_timefreq',     true  ; ...
    'PlotExponentByStage',  'exponent_by_stage',      true  ; ...
    'PlotSlopesTimecourse', 'slopes_timecourse',      true  };

suffix = '';
for k = 1:size(order, 1)
    name = order{k, 1};
    if isfield(plotflags, name) && plotflags.(name) && order{k, 3}
        suffix = order{k, 2};
    end
end
end
