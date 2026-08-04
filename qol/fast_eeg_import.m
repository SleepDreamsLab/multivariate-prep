function [EEG] = fast_eeg_import(filtFile)

    [p, f, ext] = fileparts(filtFile);
    if strcmp(ext, '.set')
        EEG = pop_loadset('filename', [f ext], 'filepath', p, 'loadmode', 'info');
        fid = fopen(fullfile(p, EEG.datfile), 'r');
        EEG.data = fread(fid, [EEG.nbchan, Inf], 'float32=>single'); fclose(fid);
        EEG = eeg_checkset(EEG);
    else
        EEG = eeg_import(filtFile);
    end

end