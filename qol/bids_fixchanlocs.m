function EEG = bids_fixchanlocs(EEG)
% QOL.BIDS_FIXCHANLOCS  Fix channel coordinates after bids_importchanlocs.
%
%   EEG = qol.bids_fixchanlocs(EEG)
%
%   bids_importchanlocs assigns X/Y/Z directly from the electrodes.tsv and
%   immediately calls eeg_checkchanlocs, which populates sph_theta from
%   X/Y/Z via cart2sph. EEGLAB/topoplot expects sph_theta=90° at the nose,
%   meaning the nose must lie along +Y. EGI electrodes.tsv uses X=anterior,
%   Y=right, so X and Y must be swapped.
%
%   Because sph_theta is already set when this function runs, a naive call
%   to eeg_checkchanlocs would reuse it and ignore the corrected X/Y. The
%   spherical fields are therefore cleared first to force recomputation.

for iCh = 1:numel(EEG.chanlocs)
    xTmp                    = EEG.chanlocs(iCh).X;
    EEG.chanlocs(iCh).X    = EEG.chanlocs(iCh).Y;   % EGI anterior → EEGLAB right? no: EGI Y=right → EEGLAB X=right
    EEG.chanlocs(iCh).Y    = xTmp;                   % EGI X=anterior → EEGLAB Y=anterior
    % Clear derived spherical fields so eeg_checkchanlocs recomputes from X/Y/Z
    EEG.chanlocs(iCh).theta      = [];
    EEG.chanlocs(iCh).radius     = [];
    EEG.chanlocs(iCh).sph_theta  = [];
    EEG.chanlocs(iCh).sph_phi    = [];
    EEG.chanlocs(iCh).sph_radius = [];
end
EEG = eeg_checkchanlocs(EEG);
end
