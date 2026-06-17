function BIDS = bidswizard(studies, base_path, rawdata_path)
arguments
    studies         = {'data-driv', 'data-elpd', 'data-ercp', 'data-ssmd', 'data-vici', 'data-wrap'}
    base_path       = 'T:\data\nin\'
    rawdata_path    = 'rawdata';
end
%   studies   : cell array of study folder names
%   base_path : string with path containing all study folders

    BIDS = cell(length(studies), 1);

    for istudy = 1:length(studies)
        current_study = studies{istudy};
        pathraw       = fullfile(base_path, current_study, rawdata_path);
        fprintf('Processing folder: %s\n', pathraw);

        tic
        % Create BIDS layout without schema validation
        BIDS{istudy} = bids.layout(pathraw, 'use_schema', 0);     
        elapsed = toc;

        fprintf('Folder %s processed in %.2f seconds.\n', current_study, elapsed);    
    end

    fprintf('Done!\n');
end
