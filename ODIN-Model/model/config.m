function cfg = config(file)
% config: parameter struct read from the Inputs sheet of assumptions.xlsx.
% file: workbook to read, assumptions.xlsx in the repository root when omitted
% cfg: struct holding all constants, orbit elements, switches, spacecraft properties and solver settings

if nargin < 1
    file = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'assumptions.xlsx');
end

C = readcell(file, 'Sheet', 'Inputs'); % C: cells of the Inputs sheet

[r0, cF] = find(cellfun(@(c) (ischar(c) || isstring(c)) && strcmp(c, 'config.m field'), C), 1);
% r0: header row of the input table
% cF: column holding the field names
if isempty(r0)
    error('config:noHeader', 'no column headed ''config.m field'' in the Inputs sheet of %s', file);
end
cV = find(cellfun(@(c) (ischar(c) || isstring(c)) && strcmp(c, 'Value'), C(r0, :)), 1); % cV: column holding the values

cfg = struct();
for k = r0+1 : size(C, 1)
    name = C{k, cF}; % name: field name in the sheet, such as cfg.drag.mass
    if ~(ischar(name) || isstring(name)) || ~startsWith(name, 'cfg.')
        continue
    end
    name = char(name);

    v = C{k, cV}; % v: value in the sheet
    if isa(v, 'missing') || isempty(v)
        error('config:noValue', '%s has no value in %s', name, file);
    end
    if isa(v, 'datetime')
        v = datevec(v);
    elseif ischar(v) || isstring(v)
        v = char(v);
        if any(strcmpi(v, {'true', 'false'}))
            v = strcmpi(v, 'true');
        elseif strcmp(name, 'cfg.epoch')
            v = datevec(v, 'yyyy-mm-dd HH:MM:SS');
        elseif ~isnan(str2double(v))
            v = str2double(v);
        end
    end

    parts = strsplit(name, '.'); % parts: cfg followed by the nested field names
    cfg = setfield(cfg, parts{2:end}, v);
end

end
