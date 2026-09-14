function cfg = defaultArtifacts()
% DEFAULTARTIFACTS  Default artifact configuration injected by bidsfun_gedai_benchmark.
%
%   cfg = bench.defaultArtifacts()
%
%   One field per artifact type (eog, emg, ecg). Numeric parameters are a fixed value or a
%   [min max] range sampled uniformly per event. Artifacts are generated OUTSIDE the cortical
%   leadfield on purpose: GEDAI keeps whatever the brain leadfield explains, so an artifact
%   projected through cortical sources would be kept by construction and the benchmark
%   would report under-cleaning that says nothing about real artifacts.
%
%   cfg.eog - corneo-retinal dipoles at both eyes, homogeneous-sphere head (bench.buildForward)
%     stages       stage digits to contaminate                        [-2 -3 0 1]
%     eyeoffset    cm from the nasion to each eye centre, in the fiducial head frame:
%                  [posterior(-)/anterior(+), lateral (+/-), inferior(-)/superior(+)]
%                                                                     [-2.5 3.2 -1.5]
%     blink        vertical rotation transient (Bell's phenomenon; the lid is not modelled
%                  separately). density /min, amplitude uV at the peak channel, duration s
%     saccade      step to a gaze angle, hold, return. density /min, angle deg, hold s,
%                  transition s, uvperdeg uV per degree at the peak channel,
%                  horizontalprob probability of a purely horizontal saccade
%     sem          slow eye movement, sinusoidal horizontal drift. density /min, angle deg,
%                  duration s, freq Hz, uvperdeg
%     psdbands     Hz, [nBands x 2] bands scored for residual power    [0.5 4]
%
%   cfg.emg - scalp-muscle bursts with local volume conduction, electrode space
%     stages, density /min, duration s
%     rms          uV RMS at the peak channel
%     band         Hz, generating band                                 [20 100]
%     sites        muscle sites; each burst picks one (temporalis, frontalis, occipital/neck,
%                  masseter)
%     sitepositions  cm, one [x y z] per site in the fiducial head frame (origin between
%                  the preauricular points, +x nasion, +y left, +z up); a burst is centred
%                  on the data channel nearest to it
%     bilateralprob  probability a temporalis/masseter burst is bilateral (e.g. a clench)
%     spread       m, Gaussian SD of spatial spread around the site electrode
%     commonfraction  variance share of a component common to all involved channels; the
%                  rest is independent per channel, so EMG is locally high-rank
%     psdbands     Hz; 45-55 Hz is skipped, where line-noise removal leaves its own notch
%
%   cfg.ecg - far-field cardiac dipole, infinite homogeneous medium (bench.buildForward)
%     stages       continuous in every epoch of these stages
%     heartposition  cm, heart centre in the fiducial head frame      [5 3 -32]
%     heartrate    beats per minute (one value per recording)          [50 70]
%     rsa          fractional RR modulation by respiration (0.25 Hz)   0.05
%     amplitude    uV, R-peak at the peak channel                      [5 20]
%     axis         electrical axis in the head frame (+x anterior, +y left, +z up)
%     axisjitter   deg, random tilt of that axis per recording
%     psdbands     Hz                                                   [1 20]
%
%   The vectorcardiogram sums Gaussian waves in cardiac phase per orthogonal lead, the
%   McSharry et al. (2003, IEEE TBME 50:289) wave model extended to three leads as in
%   Sameni et al. (2007, EURASIP J Adv Signal Process, 43407). The lead-1 wave parameters
%   are McSharry's; leads 2-3 are illustrative, chosen only to give the loop a realistic
%   rotation.

N3 = -3; N2 = -2; REM = 0; Wake = 1;
allStages = [N2 N3 REM Wake];

cfg.eog.stages    = allStages;
cfg.eog.eyeoffset = [-2.5 3.2 -1.5];
cfg.eog.blink     = struct('density', 3, 'amplitude', [100 250], 'duration', [0.25 0.45]);
cfg.eog.saccade   = struct('density', 4, 'angle', [5 25], 'hold', [0.3 1.5], ...
    'transition', [0.03 0.08], 'uvperdeg', [4 8], 'horizontalprob', 0.7);
cfg.eog.sem       = struct('density', 0.5, 'angle', [10 30], 'duration', [4 12], ...
    'freq', [0.1 0.4], 'uvperdeg', [4 8]);
cfg.eog.psdbands  = [0.5 4];

cfg.emg.stages         = allStages;
cfg.emg.density        = 2;
cfg.emg.duration       = [0.5 4];
cfg.emg.rms            = [10 40];
cfg.emg.band           = [20 100];
cfg.emg.sites          = {'temporalL', 'temporalR', 'frontalis', 'occipital', 'masseterL', 'masseterR'};
cfg.emg.sitepositions  = struct('temporalL', [1.5 7 3.5], 'temporalR', [1.5 -7 3.5], ...
    'frontalis', [9 0 3.5], 'occipital', [-9 0 -1], 'masseterL', [4.5 5.5 -4.5], 'masseterR', [4.5 -5.5 -4.5]);
cfg.emg.bilateralprob  = 0.3;
cfg.emg.spread         = 0.03;
cfg.emg.commonfraction = [0.2 0.5];
cfg.emg.psdbands       = [20 45; 55 95];

cfg.ecg.stages        = allStages;
cfg.ecg.heartposition = [5 3 -32];
cfg.ecg.heartrate     = [50 70];
cfg.ecg.rsa           = 0.05;
cfg.ecg.amplitude     = [5 20];
cfg.ecg.axis          = [0.2 0.5 -0.85];
cfg.ecg.axisjitter    = 20;
cfg.ecg.psdbands      = [1 20];
end
