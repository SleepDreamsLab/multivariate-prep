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
fprintf(fid, '%s', jsonencode(sidecar, 'PrettyPrint', true));
fclose(fid);
fprintf('[JSON] %s\n', jsonFile);
end
