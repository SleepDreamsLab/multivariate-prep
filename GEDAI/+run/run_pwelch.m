function [EpochPower, Frequencies, Times_h, KeepTime] = run_pwelch(EEG, EpochLength, WindowLength, Overlap, KeepTime)
arguments
    EEG
    EpochLength     = 30
    WindowLength    = 4
    Overlap         = .5
    KeepTime        = []
end

%%% Running pwelch
fprintf('Computing power ...\n')
D = tic;
[EpochPower, Frequencies] = oscip.compute_power_on_epochs(EEG.data, ...
    EEG.srate, EpochLength, WindowLength, Overlap);
KeepTime.Pwelch = toc(D);

%%% Save space
EpochPower  = single(EpochPower);
Frequencies = single(Frequencies);

%%% Times vector
% nEpochs = floor(EEG.pnts / EEG.srate / EpochLength);
nEpochs = size(EpochPower, 2);
Times_h = [1:nEpochs] .* EpochLength/60/60;  

%%% Report time
fprintf('Computing power took %.0fs!\n', KeepTime.Pwelch)

end