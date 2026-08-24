function sidecarjson(KeepTime, jsonFile, extra)
% WRITE_TIMING_JSON  Write a JSON sidecar with KeepTime durations converted to minutes.
%
% USAGE:
%   sidecarjson(KeepTime, jsonFile)
%   sidecarjson(KeepTime, jsonFile, struct('GEDAIParameters', r))
%
% INPUTS:
%   KeepTime  — struct with timing fields in seconds
%   jsonFile  — full path to the .json file to write
%   extra     — optional struct; fields are merged as top-level JSON keys (default empty)

arguments
    KeepTime  struct
    jsonFile  (1,1) string
    extra     struct = struct()
end

fields = fieldnames(KeepTime);
durMin = struct();
for k = 1:numel(fields)
    durMin.(fields{k}) = round(KeepTime.(fields{k}) / 60, 2);
end
sidecar = extra;
sidecar.ProcessingDurationsMinutes = durMin;
sidecar.GeneratedDate = datestr(now, 'yyyy-mm-ddTHH:MM:SS');
fid = fopen(jsonFile, 'w');
fprintf(fid, '%s', jsonencode(finitise(sidecar), 'PrettyPrint', true));
fclose(fid);
fprintf('[JSON] %s\n', jsonFile);
end

% -------------------------------------------------------------------------
function s = finitise(s)
% JSON has no Inf or NaN, and jsonencode writes both as null - so a threshold
% deliberately set to Inf to disable a criterion comes back looking like a setting
% nobody filled in. Written as the strings "Inf" / "NaN" instead, which say what was
% meant and still parse.
    if isstruct(s)
        for f = fieldnames(s)'
            for k = 1:numel(s)
                s(k).(f{1}) = finitise(s(k).(f{1}));
            end
        end
    elseif isnumeric(s) && ~isempty(s) && ~all(isfinite(s(:)))
        if isscalar(s)
            s = sprintf('%g', s);
        else
            c = num2cell(s);
            bad = ~isfinite(s);
            c(bad) = cellfun(@(v) sprintf('%g', v), c(bad), 'uni', 0);
            s = c;
        end
    end
end
