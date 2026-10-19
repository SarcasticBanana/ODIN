function ecl = eclipseTimes(out, cfg)
% eclipseTimes: shadow boundary times found by root finding between output samples.
% out: state history struct from propagate.m
% cfg: parameter struct from config.m
% ecl: struct of umbra and penumbra entry and exit times [s] and durations [s]

[gu, gp] = margins(out.r, out.rsun, cfg);
% gu: 1 x N umbra margin, negative inside the umbra [rad]
% gp: 1 x N penumbra margin, negative inside the penumbra [rad]

ecl.umbraIn = refine(gu, out, cfg, 1, -1); % ecl.umbraIn: 1 x K umbra entry times [s]
ecl.umbraOut = refine(gu, out, cfg, 1, +1); % ecl.umbraOut: 1 x K umbra exit times [s]
ecl.penumbraIn = refine(gp, out, cfg, 2, -1); % ecl.penumbraIn: 1 x K penumbra entry times [s]
ecl.penumbraOut = refine(gp, out, cfg, 2, +1); % ecl.penumbraOut: 1 x K penumbra exit times [s]

ecl.umbraDur = pairDurations(ecl.umbraIn, ecl.umbraOut); % ecl.umbraDur: 1 x K umbra durations [s]
ecl.penumbraDur = pairDurations(ecl.penumbraIn, ecl.penumbraOut); % ecl.penumbraDur: 1 x K penumbra durations [s]

end


function [gu, gp] = margins(r, rsun, cfg)
% margins: signed angular distance from the umbra and penumbra boundaries.
% r: 3 x N satellite position in ECI [km]
% rsun: 3 x N Sun position in ECI [km]
% cfg: parameter struct from config.m
% gu: 1 x N umbra margin, zero on the umbra cone [rad]
% gp: 1 x N penumbra margin, zero on the penumbra cone [rad]

[sep, aSun, aEarth] = shadowGeometry(r, rsun, cfg);

% Umbra boundary where the Earth disk just covers the Sun disk:
gu = sep - (aEarth - aSun);

% Penumbra boundary where the two disks just touch:
gp = sep - (aEarth + aSun);

end


function tz = refine(g, out, cfg, kind, dir)
% refine: brackets each sign change of the margin and solves it to solver tolerance.
% g: 1 x N margin sampled on the output grid [rad]
% out: state history struct from propagate.m
% cfg: parameter struct from config.m
% kind: 1 for the umbra boundary, 2 for the penumbra boundary
% dir: -1 for crossings into shadow, +1 for crossings into sunlight
% tz: 1 x K refined crossing times [s]

k = find(g(1:end-1).*g(2:end) < 0 & sign(g(2:end) - g(1:end-1)) == dir);
% k: 1 x K indices of the sample preceding each crossing

tz = zeros(1, numel(k));

for j = 1:numel(k)
    w = max(1, k(j)-2) : min(numel(out.t), k(j)+3); % w: local window of samples used for interpolation
    tz(j) = fzero(@(tt) marginAt(tt, out.t(w), out.r(:,w), cfg, kind), ...
                  [out.t(k(j)) out.t(k(j)+1)]);
end

end


function m = marginAt(tt, tw, rw, cfg, kind)
% marginAt: margin at an arbitrary time, position splined and Sun vector evaluated exactly.
% tt: time since sim start [s]
% tw: 1 x W window of sample times [s]
% rw: 3 x W window of sample positions in ECI [km]
% cfg: parameter struct from config.m
% kind: 1 for the umbra boundary, 2 for the penumbra boundary
% m: margin at tt [rad]

r = interp1(tw, rw.', tt, 'spline').'; % r: interpolated position in ECI [km]
rsun = sunVector(tt, cfg); % rsun: Sun position in ECI [km]

[gu, gp] = margins(r, rsun, cfg);

if kind == 1
    m = gu;
else
    m = gp;
end

end


function d = pairDurations(tIn, tOut)
% pairDurations: time from each entry to the exit that follows it.
% tIn: 1 x K entry times [s]
% tOut: 1 x K exit times [s]
% d: 1 x K durations [s]

d = zeros(1, numel(tIn));

for j = 1:numel(tIn)
    n = find(tOut > tIn(j), 1); % n: index of the first exit after this entry
    if isempty(n)
        d(j) = NaN;
    else
        d(j) = tOut(n) - tIn(j);
    end
end

end
