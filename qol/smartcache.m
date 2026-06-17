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
        fileID = [fileID '.mat'];
    else
        fileID = [fileID ext];
    end

    % skipMask: '' entries are collected from func but not saved or returned
    skipMask  = cellfun(@isempty, variablenames);
    saveNames = variablenames(~skipMask);
    nCollect  = numel(variablenames);   % positional outputs to request from func
    nReturn   = numel(saveNames);       % outputs to save and expose to caller

    %%% Load from cache
    if ~refresh && isfile(filename)
        D = tic;
        fprintf('Loading %s ...\n', filename)
        S = load(filename);
        fprintf('Done!\n'); fprintf('Loading took %.0fs!\n', toc(D))
        S = orderfields(S, saveNames);
        varargout = struct2cell(S);
        return
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

    %%% Save cache
    fprintf('Saving %s in %s ...\n', fileID, filepath)
    if whos('out').bytes > 2e9
        save(filename, '-fromstruct', out, '-v7.3')
    else
        save(filename, '-fromstruct', out)
    end
    fprintf('Done!\n')
end
