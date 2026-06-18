function [Results, KeepTime] = run_fooof(SmoothPower, Frequencies, KeepTime, varargin)

%%% Default variables
if nargin < 3 || isempty(KeepTime)
    KeepTime = [];
end

%%% Default parameters
opts.FrequencyRange      = [2 45];
% opts.MeanSmoothSpan      = 2;
% opts.MedianSmoothSpan    = 3;
opts.peak_width_limits   = [0.5 20];
opts.aperiodic_mode      = 'fixed';

% Parse name-value pairs (varargin)
for k = 1:2:numel(varargin)
    opts.(varargin{k}) = varargin{k+1};
end

%%% ***********************************************************************



% % smooth signal so FOOOF works more smoothly
% A = tic;
% SmoothPower = oscip.smooth_spectrum_median(EpochPower, Frequencies, opts.MedianSmoothSpan); 
% SmoothPower = oscip.smooth_spectrum(SmoothPower, Frequencies, opts.MeanSmoothSpan); 
% KeepTime.Smoothing = toc(A);

% For Specparam
AdditionalParameters = struct();
AdditionalParameters.peak_width_limits  = opts.peak_width_limits;
AdditionalParameters.aperiodic_mode     = opts.aperiodic_mode;

% run specparam
B = tic;
% [Exponents, Offsets, FrequenciesPeriodic, PeriodicPeaks, PeriodicPower, Errors, RSquared, Knees] = ...
%     oscip.fit_fooof_multidimentional_matlab(SmoothPower, Frequencies, opts.FrequencyRange, AdditionalParameters);
[Exponents, Offsets, FrequenciesPeriodic, PeriodicPeaks, PeriodicPower, Errors, RSquared] = ...
    oscip.fit_fooof_multidimentional_matlab(SmoothPower, Frequencies, opts.FrequencyRange, AdditionalParameters);
KeepTime.FOOOF = toc(B);

% turn everything to singles
% Results.EpochPower      = single(EpochPower);
% Results.SmoothPower     = single(SmoothPower);
Results.Exponents       = single(Exponents);
Results.Offsets         = single(Offsets);
Results.FrequenciesPeriodic = single(FrequenciesPeriodic);
Results.PeriodicPeaks   = single(PeriodicPeaks);
Results.PeriodicPower   = single(PeriodicPower);
Results.Errors          = single(Errors);
Results.RSquared        = single(RSquared);
% Results.Knees           = single(Knees);
Results.KeepTime        = KeepTime;

end
