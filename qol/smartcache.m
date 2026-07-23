function varargout = smartcache(func, filename, refresh, variablenames, varargin)
% SMARTCACHE  Run func() and cache its outputs to a .mat file.
%
%   Use '' as a placeholder in variablenames to collect a positional output
%   from func without saving it or returning it.  E.g.:
%     [a, c, d] = smartcache(f, file, refresh, {'a','','c','d'})
%   collects 4 outputs from f but only saves/returns a, c, d.

    %%% Check filename extension
    [filepath, fileID, ext] = fileparts(filename);
    if isempty(ext)
        ext   = '.mat';
        fileID = [fileID ext];
    else
        fileID = [fileID ext];
    end
    isBrainvision = strcmpi(ext, '.vhdr');

%     %%% Truncate filename if full path exceeds Windows MAX_PATH (260 chars)
%     MAX_PATH = 260;
%     if numel(filename) > MAX_PATH
%         allowedStem = MAX_PATH - numel(filepath) - numel(ext) - 1; % -1 for filesep
%         if allowedStem < 1
%             error('smartcache:pathTooLong', ...
%                 'Directory path alone (%d chars) exceeds MAX_PATH=%d.', ...
%                 numel(filepath), MAX_PATH)
%         end
%         [~, stemOnly] = fileparts(fileID);   % fileID already has ext; strip it
%         newStem  = stemOnly(1:allowedStem);
%         fileID   = [newStem ext];
%         filename = fullfile(filepath, fileID);
%         warning('smartcache:pathTruncated', ...
%             'Path exceeded %d chars; filename stem truncated to:\n  %s', MAX_PATH, filename)
%     end

    % skipMask: '' entries are collected from func but not saved or returned
    skipMask  = cellfun(@isempty, variablenames);
    saveNames = variablenames(~skipMask);
    nCollect  = numel(variablenames);   % positional outputs to request from func
    nReturn   = numel(saveNames);       % outputs to save and expose to caller

    %%% Load from cache
    cacheOk = ~refresh && isfile(filename) && getfield(dir(filename), 'bytes') > 0;
    if cacheOk
        D = tic;
        fprintf('Loading %s ...\n', filename)
        try
            if isBrainvision
                varargout = {eeg_import(filename)};
            else
                S = load(filename);
                S = orderfields(S, saveNames);
                varargout = struct2cell(S);
            end
            fprintf('Done!\n'); fprintf('Loading took %.0fs!\n', toc(D))
            return
        catch ME
            warning('smartcache:loadFailed', ...
                'Cache file appears corrupt, re-running: %s', ME.message)
        end
    end

    %%% Run function — always request nCollect outputs to reach every position
    tmp = cell(1, nCollect);
    [tmp{1:nCollect}] = func();

    %%% Pack non-skipped outputs into struct and varargout
    out     = struct();
    iReturn = 0;
    varargout = cell(1, nReturn);
    for k = 1:nCollect
        if ~skipMask(k)
            iReturn = iReturn + 1;
            out.(variablenames{k}) = tmp{k};
            varargout{iReturn}     = tmp{k};
        end
    end

    %%% Create folder if needed
    if ~isempty(filepath) && ~exist(filepath, 'dir')
        mkdir(filepath)
    end

    %%% Save cache  (delete any stale/0-byte file first so save cannot see a "corrupt" file)
    fprintf('Saving %s in %s ...\n', fileID, filepath)
    if isfile(filename)
        delete(filename)
    end
    if isBrainvision
        pop_writebva(varargout{1}, filename, 'DataOrientation', 'MULTIPLEXED');
    elseif whos('out').bytes > 2e9
        save(filename, '-fromstruct', out, '-v7.3')
    else
        save(filename, '-fromstruct', out)
    end
    fprintf('Done!\n')
end
