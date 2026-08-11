function [EEG, chanmap] = chans1020(EEG, select, opts)
arguments
    EEG
    select = false;
    opts.net {mustBeMember(opts.net, {'EGI256', 'EGI128'})} = 'EGI256';
    opts.add_eog = false
    opts.chanprefix char = ''   % '' = chanmap numbers are chanlocs indices; e.g. 'E' = look up channels by label 'E<n>'
end

%%% EGI 256 -> 10-20 mapping  (Cz omitted – absent in standard EGI table)
switch opts.net
    case 'EGI256'
        chanmap = struct( ...
            'Fp1',  37,  ...
            'Fpz',  26,  ...
            'Fp2',  18,  ...
            'F7',   47,  ...
            'F3',   36,  ...
            'Fz',   21,  ...
            'F4',  224,  ...
            'F8',    2,  ...
            'T7',   69,  ...
            'C3',   59,  ...
            'C4',  183,  ...
            'T8',  202,  ...
            'P7',   96,  ...
            'P3',   87,  ...
            'Pz',  101,  ...
            'P4',  153,  ...
            'P8',  170,  ...
            'O1',  116,  ...
            'O2',  150   ...
        );
        eogchans = [54 248 1 230];

    case 'EGI128'
        % 124 channels
        chanmap = struct( ...
            'Fp1',  22,   'Fp2',  9, 'F7',   33,  'F3',   24,  'Fz',   11,  'F4', 120, 'F8',   118, ...
             'C3',   36, 'C4', 102,  'P3',   51,  'Pz',  60,  'P4', 90, ...
            'O1',  68,  'O2',  81);
            eogchans = [128-4 32 1 125-4];

end

% Structure content
labels1020  = fieldnames(chanmap);
enum1020    = struct2array(chanmap);    % original channel numbers (as used by '<prefix><n>' labels)

% Resolve 10-20 channels to their current position in EEG.chanlocs.
% chanprefix == '': chanmap numbers are assumed to still be valid chanlocs
% indices (i.e. no channels removed upstream). Otherwise, channels are
% looked up by '<chanprefix><n>' label (e.g. 'E22'), so it still works if
% bad channels were previously dropped and indices shifted.
if isempty(opts.chanprefix)
    idx1020 = enum1020;
else
    idx1020 = chanNumToIndex(EEG.chanlocs, enum1020, opts.chanprefix);
end

% Relabel (channels removed upstream and not found simply stay unlabelled)
for iCh = 1:numel(idx1020)
    if ~isnan(idx1020(iCh))
        EEG.chanlocs(idx1020(iCh)).labels = labels1020{iCh};
    end
end

% Add EOG
% Bipolar EOG: EOG1 = ch54 - ch1,  EOG2 = ch230 - ch248
if opts.add_eog
    if isempty(opts.chanprefix)
        eogData = EEG.data(eogchans, :);
    else
        eogData = resolveEogData(EEG, eogchans, opts.chanprefix);
    end
    EEG.data(end+1,:) = eogData(1,:) - eogData(3,:);
    EEG.chanlocs(end+1).labels = 'EOG1';
    EEG.data(end+1,:) = eogData(4,:) - eogData(2,:);
    EEG.chanlocs(end+1).labels = 'EOG2';
    EEG.nbchan = size(EEG.data, 1);
    idx1020 = [idx1020, EEG.nbchan-1, EEG.nbchan];
end

% Select 10-20 electrodes (channels not found upstream, i.e. NaN, are dropped)
if select
    EEG = pop_select(EEG, 'channel', idx1020(~isnan(idx1020)));
end

% EOG1/EOG2 have no '<prefix><n>' label of their own (they're synthetic,
% not part of the net), so record their actual channel index directly.
if opts.add_eog
    chanmap.EOG1 = EEG.nbchan - 1;
    chanmap.EOG2 = EEG.nbchan;
end

end % chans1020

% -------------------------------------------------------------------------
function idx = chanNumToIndex(chanlocs, chanNum, chanprefix)
% Map original '<chanprefix><n>' channel numbers to their current index in
% chanlocs by label. Returns NaN for numbers that are no longer present
% (e.g. removed as bad channels upstream).
    wantLabels = arrayfun(@(n) sprintf('%s%d', chanprefix, n), chanNum, 'uni', 0);
    curLabels  = {chanlocs.labels};
    idx = nan(1, numel(chanNum));
    for iCh = 1:numel(chanNum)
        pos = find(strcmp(curLabels, wantLabels{iCh}), 1);
        if ~isempty(pos), idx(iCh) = pos; end
    end
end

% -------------------------------------------------------------------------
function eogData = resolveEogData(EEG, eogchans, chanprefix)
% Fetch the 4 EOG anchor channels ('<chanprefix><n>' labels) by label. Any
% that were removed upstream are interpolated from the surviving channels
% using EEG.urchanlocs (which retains all original channel locations) via
% spherical-spline interpolation.
    eogLabels = arrayfun(@(n) sprintf('%s%d', chanprefix, n), eogchans, 'uni', 0);
    curLabels = {EEG.chanlocs.labels};
    present   = ismember(eogLabels, curLabels);

    EEGext = EEG;
    if ~all(present)
        if ~isfield(EEG, 'urchanlocs') || isempty(EEG.urchanlocs)
            error('chans1020:noUrchanlocs', ...
                'EEG.urchanlocs is required to interpolate missing EOG anchor channel(s).');
        end
        urLabels = {EEG.urchanlocs.labels};
        [found, urIdx] = ismember(eogLabels(~present), urLabels);
        if ~all(found)
            error('chans1020:eogNotFound', ...
                'Missing EOG anchor channel(s) not found in EEG.urchanlocs either.');
        end
        EEGext = eeg_interp(EEGext, EEG.urchanlocs(urIdx(found)), 'spherical');
    end

    extLabels = {EEGext.chanlocs.labels};
    [~, pos]  = ismember(eogLabels, extLabels);
    eogData   = EEGext.data(pos, :);
end
