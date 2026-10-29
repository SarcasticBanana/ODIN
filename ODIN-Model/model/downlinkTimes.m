function [dl, con] = downlinkTimes(out, cfg, st)
% downlinkTimes: start and end times of the windows in which the satellite is above the mask of at least one ground station.
% out: state history struct from propagate.m
% cfg: parameter struct from config.m
% st: 1 x M station structs from stations.m, the stations in use in the Stations sheet when omitted
% dl: K x 2 start and end times since sim start, one row per window, windows of different stations merged where they overlap, cut at the run ends [s]
% con: 1 x M cell of contact window structs from contactTimes.m, one per station

if nargin < 3
    st = stations();
end

con = cell(1, numel(st));
iv = zeros(0, 2); % iv: windows of every station, cut at the run ends [s]
for g = 1:numel(st)
    con{g} = contactTimes(out, cfg, st(g));
    iv = [iv; max(con{g}.aos(:), out.t(1)), min(con{g}.los(:), out.t(end))];
end
iv = sortrows(iv);

% Windows that overlap merged into one:
dl = zeros(0, 2);
for k = 1:size(iv, 1)
    if ~isempty(dl) && iv(k, 1) <= dl(end, 2)
        dl(end, 2) = max(dl(end, 2), iv(k, 2));
    else
        dl(end+1, :) = iv(k, :);
    end
end

end
