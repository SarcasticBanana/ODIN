function con = contactTimes(out, cfg, st)
% contactTimes: windows in which a ground station sees the satellite above its elevation mask, edges found by root finding between output samples.
% out: state history struct from propagate.m
% cfg: parameter struct from config.m
% st: station struct from stations.m
% con: struct of acquisition and loss times [s], highest elevation [deg] and its time [s] per window, NaN where a window runs past either end of the run

N = numel(out.t); % N: number of output samples
g = margin(out.r, out.t, st, cfg); % g: 1 x N elevation above the mask at the output samples [deg]

kUp = find(g(1:end-1) < 0 & g(2:end) >= 0); % kUp: samples preceding each rise through the mask
kDn = find(g(1:end-1) >= 0 & g(2:end) < 0); % kDn: samples preceding each fall through the mask
gx = [-Inf, g, -Inf]; % gx: margins padded so the first and last samples can be maxima [deg]
kPk = find(gx(2:end-1) < 0 & gx(2:end-1) >= gx(1:end-2) & gx(2:end-1) > gx(3:end)); % kPk: sampled maxima below the mask, where a pass shorter than a sample can sit

tUp = zeros(1, numel(kUp)); % tUp: rise times [s]
for j = 1:numel(kUp)
    w = max(1, kUp(j) - 2) : min(N, kUp(j) + 3); % w: samples splined around the crossing
    tUp(j) = fzero(@(x) marginAt(x, out, w, st, cfg), out.t(kUp(j) + [0 1]));
end

tDn = zeros(1, numel(kDn)); % tDn: fall times [s]
for j = 1:numel(kDn)
    w = max(1, kDn(j) - 2) : min(N, kDn(j) + 3);
    tDn(j) = fzero(@(x) marginAt(x, out, w, st, cfg), out.t(kDn(j) + [0 1]));
end

for k = kPk
    w = max(1, k - 3) : min(N, k + 3);
    tb = out.t([max(1, k - 1), min(N, k + 1)]); % tb: neighbouring sample times bracketing the maximum [s]
    [tPk, gPk] = fminbnd(@(x) -marginAt(x, out, w, st, cfg), tb(1), tb(2));
    % tPk: time of the highest elevation between the neighbouring samples [s]
    % gPk: minus the elevation above the mask there [deg]
    if -gPk > 0
        tUp(end+1) = fzero(@(x) marginAt(x, out, w, st, cfg), [tb(1), tPk]);
        tDn(end+1) = fzero(@(x) marginAt(x, out, w, st, cfg), [tPk, tb(2)]);
    end
end

ev = sortrows([tUp(:), ones(numel(tUp), 1); tDn(:), -ones(numel(tDn), 1)]); % ev: event times [s] with +1 for a rise and -1 for a fall
aos = zeros(1, 0); % aos: acquisition times [s]
los = zeros(1, 0); % los: loss times [s]
if g(1) >= 0
    aos(end+1) = NaN;
end
for j = 1:size(ev, 1)
    if ev(j, 2) > 0
        aos(end+1) = ev(j, 1);
    else
        los(end+1) = ev(j, 1);
    end
end
if numel(los) < numel(aos)
    los(end+1) = NaN;
end

K = numel(aos); % K: number of windows
con.aos = aos; % con.aos: 1 x K acquisition times since sim start, NaN when the run starts in view [s]
con.los = los; % con.los: 1 x K loss times since sim start, NaN when the run ends in view [s]
con.elMax = zeros(1, K); % con.elMax: 1 x K highest elevation of each window [deg]
con.tMax = zeros(1, K); % con.tMax: 1 x K time of the highest elevation [s]

for j = 1:K
    a = max(aos(j), out.t(1)); % a: window start inside the run [s]
    b = min(los(j), out.t(end)); % b: window end inside the run [s]
    if isnan(a)
        a = out.t(1);
    end
    if isnan(b)
        b = out.t(end);
    end
    w = max(1, find(out.t <= a, 1, 'last') - 3) : min(N, find(out.t >= b, 1) + 3);
    [con.tMax(j), gMax] = fminbnd(@(x) -marginAt(x, out, w, st, cfg), a, b); % gMax: minus the highest elevation above the mask [deg]
    con.elMax(j) = st.mask - gMax;
end

end


function g = margin(r, t, st, cfg)
% margin: elevation of the satellite above the station's mask.
% r: 3 x N satellite position in the model frame [km]
% t: 1 x N times since sim start [s]
% st: station struct from stations.m
% cfg: parameter struct from config.m
% g: 1 x N elevation above the mask [deg]

[rgs, up] = stationVector(t, st, cfg);
% rgs: 3 x N station position [km]
% up: 3 x N local vertical at the station

rho = r - rgs; % rho: 3 x N line of sight from the station to the satellite [km]

% Elevation, the angle between the line of sight and the local horizontal plane:
g = asind(sum(rho.*up, 1)./vecnorm(rho)) - st.mask;

end


function g = marginAt(x, out, w, st, cfg)
% marginAt: elevation above the mask at an arbitrary time, position splined across a window of samples.
% x: time since sim start [s]
% out: state history struct from propagate.m
% w: indices of the samples splined
% st: station struct from stations.m
% cfg: parameter struct from config.m
% g: elevation above the mask at x [deg]

r = interp1(out.t(w), out.r(:, w).', x, 'spline').'; % r: satellite position at x [km]
g = margin(r, x, st, cfg);

end
