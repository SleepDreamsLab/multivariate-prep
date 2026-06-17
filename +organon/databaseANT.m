%%% General parameters
Parameters.EpochLength                              = 4;       % in seconds
Parameters.StageMap                                 = dictionary( ...
        'W',    0, ...
        'REM',  5, ...
        'N1',   1, ...
        'N2',   2, ...
        'N3',   3 ...
        );

% %%% Channels parameters
% Parameters.Channels.Mastoids    = [94 190];
% Parameters.Channels.NotEEG      = [258:264]; % 133 145 165 174 
% Parameters.Channels.Face        = [67 73 82 91 216:219 225:256];
% Parameters.Channels.EEG         = 1:264;
% Parameters.Channels.EEG([ ...
%     Parameters.Channels.NotEEG, ...
%     Parameters.Channels.Face ...
%     ]) = [];
% 
% % Channels that really should have good signal quality; good for averaging
% Parameters.Channels.Neck        = [92 93 102:104 111 112 120 121 133 134 145 ... 
%     146 156 165 166 174 175 187:189 199:201 208 209];
% Parameters.Channels.Edge        = [105 114 123 136 148 158 168 177];
% Parameters.Channels.CoreEEG     = 1:264;
% Parameters.Channels.CoreEEG([ ...
%     Parameters.Channels.NotEEG, ...
%     Parameters.Channels.Face, ...
%     Parameters.Channels.Edge, ...
%     Parameters.Channels.Neck ...
%     ]) = [];

% Channels for topography
Parameters.Channels.TopoChannels = 1:264;
Parameters.Channels.TopoChannels([Parameters.Channels.NotEEG]) = [];

% Channel parameters actually used
Parameters.Channels.SelectedChannels                = 1:64;  % channel indices
Parameters.Channels.MaxChannelsToInterpolate        = 20;       % # channels
Parameters.Channels.MaxBadMuscleChannels            = 10;       % # channels

%%% Interpolation
Parameters.Interpolation.MinGapEpochs               = 2;        % in epochs
Parameters.Interpolation.MinEpochs                  = 30;       % in epochs


%%% Pwelch parameters
Parameters.Welch.WindowLength                       = 4;        % in seconds
Parameters.Welch.Overlap                            = 0.5;      % in proportion

%%% FOOF parameters
Parameters.FOOOF.FrequencyRange                     = [2 45];   % in Hz (considered range)
Parameters.FOOOF.MeanSmoothSpan                     = 2;        % in Hz
Parameters.FOOOF.MedianSmoothSpan                   = 3;        % in Hz
Parameters.FOOOF.PeakWidthLimits                    = [.5 20];
Parameters.FOOOF.AperiodicMode                      = 'fixed';  % 'fixed', 'knee'
Parameters.FOOOF.ExponentVarianceThreshold          = 0.6;

%%% ICA parameters
Parameters.ICA.ICAMinutes                           = [2 15];   % in minutes (min & max)
Parameters.ICA.ICAArtefactTypes                     = {'Eye', 'Heart', 'Muscle'};
Parameters.ICA.MaxBadMuscleChannels                 = 10;       % channel number
Parameters.ICA.MinTimeSleepCycles                   = 10;       % in minutes; to define a sleep cycle, need at least these many minutes of NREM data

%%% Thresholds: artefact in time domain
Parameters.ThreshTime.Big.MaxAmplitude                          = 500;      % in microvolt
Parameters.ThreshTime.Big.MaxDiff                               = 100;      % in microvolt per second
Parameters.ThreshTime.Big.PaddingLength                         = 1;        % in seconds
Parameters.ThreshTime.Flat.MinVariance                          = 1.8;      % in microvolt
Parameters.ThreshTime.Flat.MovingVarianceWindowLength           = 10;       % in seconds
Parameters.ThreshTime.Detached.MinVariance                      = 2;        % in seconds
Parameters.ThreshTime.Detached.MovingVarianceWindowLength       = 300;      % in seconds
Parameters.ThreshTime.Muscle.FrequencyRange                     = [20 95];  % in Hz
Parameters.ThreshTime.Muscle.MedianMultiplierThresholds         = [20 100]; % in ?
Parameters.ThreshTime.Muscle.SmoothWindow                       = .2        % in seconds
Parameters.ThreshTime.Disconnected.CorrelationWindow            = 30;       % in seconds
Parameters.ThreshTime.Disconnected.MaxCorrelationThreshold      = .999;
Parameters.ThreshTime.Disconnected.MinCorrelationThreshold      = .3;
Parameters.ThreshTime.Disconnected.MinCorrChannels              = 5;
Parameters.ThreshTime.Correlation.CorrelationThreshold          = .3;   % theoretically from [-1 1]; recommended between 0-.5
Parameters.ThreshTime.Correlation.DifferenceThreshold           = 5;    % how many times more than the standard deviation of the neighboring channel is still acceptable
Parameters.ThreshTime.Correlation.HighPassFilter                = 0.8;  % ignore drifts and sweating artefacts
Parameters.ThreshTime.Correlation.LowPassFilter                 = 6;    % probably good from 4- 15 Hz; higher frequency activity doesn't correlate so much between neighbors
Parameters.ThreshTime.Correlation.CorrWindow                    = 4;    % in seconds. The longer the window, the less difference there is between bad channels and actually good channels
Parameters.ThreshTime.Correlation.STDWindow                     = 60*5; % to calculate moving standard deviation; uses a long window so less affected by occasional artefacts

%%% Thresholds: artefact in power domain
Parameters.ArtefactsSpecparam.Exponent.Range                           = [0.5 4];
Parameters.ArtefactsSpecparam.Exponent.OutlierMovingWindow             = 5;        % in minutes
Parameters.ArtefactsSpecparam.Exponent.OutlierMaxDeviation             = .75;
Parameters.ArtefactsSpecparam.Exponent.VarianceThreshold               = .6;
Parameters.ArtefactsSpecparam.Offset.Range                             = [-1 3.5]; % Sophia: [-1 5]
Parameters.ArtefactsSpecparam.Offset.OutlierMovingWindow               = 5;        % in minutes
Parameters.ArtefactsSpecparam.Offset.OutlierMaxDeviation               = 1.5;   
Parameters.ArtefactsSpecparam.Offset.VarianceThreshold                 = .6;
Parameters.ArtefactsSpecparam.Error.MaxError                           = .15;
Parameters.ArtefactsSpecparam.RSquared.MinRSquared                     = .95;
Parameters.ArtefactsSpecparam.ExpOffsetRelation.zResidualThreshold     = 10; % Sophia: 10;  how many z-scores away from the correlation line between exponents and offsets can a datapoint be before removed
