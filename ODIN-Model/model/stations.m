function st = stations(file)
% stations: ground stations marked for use in the Stations sheet of assumptions.xlsx.
% file: workbook to read, assumptions.xlsx in the repository root when omitted
% st: 1 x M struct array of station name, geodetic latitude [deg], longitude [deg], altitude above the ellipsoid [m] and elevation mask [deg]

if nargin < 1
    file = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'assumptions.xlsx');
end

C = readcell(file, 'Sheet', 'Stations'); % C: cells of the Stations sheet

[r0, cN] = find(cellfun(@(c) (ischar(c) || isstring(c)) && strcmp(c, 'Name'), C), 1);
% r0: header row of the station table
% cN: column holding the station names
if isempty(r0)
    error('stations:noHeader', 'no column headed ''Name'' in the Stations sheet of %s', file);
end

hdr = C(r0, :); % hdr: header row cells
col = @(s) find(cellfun(@(c) (ischar(c) || isstring(c)) && startsWith(c, s), hdr), 1); % col: column whose header starts with the given text
cLat = col('Latitude'); % cLat: latitude column
cLon = col('Longitude'); % cLon: longitude column
cAlt = col('Altitude'); % cAlt: altitude column
cMask = col('Mask'); % cMask: elevation mask column
cUse = col('Use'); % cUse: use flag column

st = struct('name', {}, 'lat', {}, 'lon', {}, 'alt', {}, 'mask', {});
for k = r0+1 : size(C, 1)
    name = C{k, cN}; % name: station name in the sheet
    if ~(ischar(name) || isstring(name)) || ~isUsed(C{k, cUse})
        continue
    end
    v = C(k, [cLat cLon cAlt cMask]); % v: latitude, longitude, altitude and mask cells of this station
    if ~all(cellfun(@(c) isnumeric(c) && isscalar(c) && isfinite(c), v))
        error('stations:noValue', 'station %s needs numeric latitude, longitude, altitude and mask in %s', char(name), file);
    end
    st(end+1) = struct('name', char(name), 'lat', v{1}, 'lon', v{2}, 'alt', v{3}, 'mask', v{4});
end

end


function tf = isUsed(c)
% isUsed: true for a ticked use cell, a logical true or the text 'true'.
% c: use cell from the Stations sheet
% tf: true when the station is in use

tf = (islogical(c) && isscalar(c) && c) || ((ischar(c) || isstring(c)) && strcmpi(c, 'true'));

end
