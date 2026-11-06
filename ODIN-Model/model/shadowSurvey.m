function sv = shadowSurvey(cfg, nDay)
% shadowSurvey: Earth shadow of every orbit over nDay days of the ode45 orbit, the longest each day, the Sun angle and LTAN each day, the node rate, nodal period and orbit-mean altitude, and the state history over the first week and over the week ending 0.7 nodal periods after the longest shadow.
% cfg: parameter struct from config.m, drag as set there, off to hold the altitude
% nDay: days the shadows are counted over from sim start [days]
% sv: struct, each field described where it is set

cfg.use.eclipse = true;

% Kepler's third law:
T = 2*pi*sqrt((cfg.Re + cfg.h)^3/cfg.mu); % T: orbital period [s]

out = propagate(cfg, (nDay + 1)*86400/T); % out: state history, one day past nDay so the orbit after the last counted shadow is there
ecl = eclipseTimes(out, cfg); % ecl: umbra and penumbra entry and exit times from eclipseTimes.m [s]

% Shadows entered within the counted days, penumbra entry to exit:
q = ecl.penumbraIn < nDay*86400 & ~isnan(ecl.penumbraDur); % q: shadows counted
sv.in = ecl.penumbraIn(q); % sv.in: 1 x K shadow entry times [s]
sv.dur = ecl.penumbraDur(q); % sv.dur: 1 x K shadow lengths, penumbra entry to exit [s]

d = floor(sv.in/86400) + 1; % d: 1 x K day of each shadow, 1 on the first day
sv.day = accumarray(d(:), sv.dur(:), [nDay 1], @max, 0).'; % sv.day: 1 x nDay longest shadow of each day, zero without shadow [s]

tD = (0:nDay-1)*86400; % tD: 1 x nDay day starts, 00:00 of each day from sim start [s]
iD = round(tD/cfg.ode.dtOut) + 1; % iD: 1 x nDay output samples at the day starts
hn = out.hvec(:, iD)./vecnorm(out.hvec(:, iD)); % hn: 3 x nDay unit orbit normal at the day starts

% Sun elevation above the orbit plane, beta = asin(h_hat . s_hat):
sv.beta = asind(sum(hn.*out.sunhat(:, iD), 1)); % sv.beta: 1 x nDay beta angle at each day start [deg]

raan = rad2deg(unwrap(atan2(out.hvec(1, iD), -out.hvec(2, iD)))); % raan: 1 x nDay right ascension of the ascending node, unwrapped [deg]
jd0 = datenum(cfg.epoch) + 1721058.5; % jd0: Julian date at sim start [days]

% Mean Sun right ascension as in tradeLighting.m, L = 280.460 + 36000.771 T deg:
raSun = 280.460 + 36000.771*(jd0 + tD/86400 - 2451545.0)/36525; % raSun: 1 x nDay mean Sun right ascension [deg]

% Local time of the ascending node, LTAN = 12 h + (RAAN - RA_sun)/(15 deg/h):
sv.ltan = mod(12 + (raan - raSun)/15, 24); % sv.ltan: 1 x nDay local time of the ascending node at each day start [h]

% Node rate, least-squares line through the daily node:
pW = polyfit(tD/86400, raan, 1); % pW: slope [deg/day] and intercept [deg]
sv.dW = pW(1); % sv.dW: node rate of the ode45 orbit [deg/day]

% Ascending node crossings over the first three days, z from below to above zero:
k3 = find(out.t <= 3*86400, 1, 'last'); % k3: last sample of the first three days
z = out.r(3, 1:k3); % z: 1 x k3 z component of position [km]
kc = find(z(1:end-1) <= 0 & z(2:end) > 0); % kc: samples before each crossing
tc = out.t(kc) - z(kc).*(out.t(kc+1) - out.t(kc))./(z(kc+1) - z(kc)); % tc: crossing times, linear between samples [s]
sv.TN = (tc(end) - tc(1))/(numel(tc) - 1); % sv.TN: nodal period, time between ascending node crossings [s]

% Orbit-mean radius over whole nodal periods, trapezoidal with the end pieces to the crossings:
rn = vecnorm(out.r(:, 1:k3)); % rn: 1 x k3 distance from the centre of the Earth [km]
i1 = kc(1) + 1; % i1: first sample after the first crossing
i2 = kc(end); % i2: last sample before the last crossing
r1 = interp1(out.t(1:k3), rn, tc(1)); % r1: radius at the first crossing [km]
r2 = interp1(out.t(1:k3), rn, tc(end)); % r2: radius at the last crossing [km]
I = trapz(out.t(i1:i2), rn(i1:i2)) + (out.t(i1) - tc(1))*(r1 + rn(i1))/2 + (tc(end) - out.t(i2))*(r2 + rn(i2))/2; % I: integral of the radius [km s]
sv.h0 = I/(tc(end) - tc(1)) - cfg.Re; % sv.h0: orbit-mean altitude above cfg.Re over the first three days [km]

% Longest shadow over the counted days:
sv.dMax = 0; % sv.dMax: longest shadow, zero without shadow [s]
sv.tMax = NaN; % sv.tMax: entry time of the longest shadow [s]
sv.betaMax = NaN; % sv.betaMax: beta angle at that entry [deg]
sv.ltanMax = NaN; % sv.ltanMax: LTAN at that entry [h]
if ~isempty(sv.dur)
    [sv.dMax, iM] = max(sv.dur); % iM: index of the longest shadow
    sv.tMax = sv.in(iM);
    kM = round(sv.tMax/cfg.ode.dtOut) + 1; % kM: output sample nearest the entry
    sv.betaMax = asind(dot(out.hvec(:, kM)/norm(out.hvec(:, kM)), out.sunhat(:, kM)));
    raM = 280.460 + 36000.771*(jd0 + out.t(kM)/86400 - 2451545.0)/36525; % raM: mean Sun right ascension at the entry [deg]
    sv.ltanMax = mod(12 + (atan2d(out.hvec(1, kM), -out.hvec(2, kM)) - raM)/15, 24);
end

% Weeks for the power runs, the worst one ending 0.7 nodal periods after the longest shadow entry so its last orbit holds that shadow:
w = 7*86400; % w: week length [s]
sv.first = slice(out, 0, w); % sv.first: state history over the first week
tE = min(sv.tMax + 0.7*sv.TN, out.t(end)); % tE: end of the worst week [s]
if isnan(tE)
    tE = w;
end
sv.worst = slice(out, max(tE - w, 0), tE); % sv.worst: state history over the worst week

end


function o = slice(out, t0, t1)
% slice: state history between two times, every history field cut to the samples inside.
% out: state history struct from propagate.m
% t0: start time [s]
% t1: end time [s]
% o: the same struct over t0 to t1

q = out.t >= t0 - 1e-6 & out.t <= t1 + 1e-6; % q: samples inside
o = out;
f = {'t', 'r', 'v', 'hvec', 'rsun', 'sunhat', 'illum'}; % f: history fields
for j = 1:numel(f)
    o.(f{j}) = out.(f{j})(:, q);
end

end
