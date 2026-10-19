function res = tradeLighting(nDays, idx)
% tradeLighting: propagates the mission orbit trade cases and finds their exact shadow intervals.
% nDays: propagation length from cfg.epoch [days]
% idx: indices of the cases to run, 1 to 6 for A to F
% res: struct array of case name, label, parameters, orbital period, run length and shadow times

if nargin < 1
    nDays = 365;
end
if nargin < 2
    idx = 1:6;
end

base = config(); % base: parameter struct shared by every case

name = {'A', 'B', 'C', 'D', 'E', 'F'}; % name: case letters from the orbit trade study
hp = [650 650 800 950 800 650]; % hp: perigee altitude above base.Re [km]
ha = [650 650 800 950 800 950]; % ha: apogee altitude above base.Re [km]
ltan = [18 10.5 18 18 18 18]; % ltan: mean local solar time of the ascending node at epoch [h]
incFix = [NaN NaN NaN NaN 90 NaN]; % incFix: fixed inclination, NaN for Sun synchronous [deg]

% Julian date from MATLAB serial date number, JD = datenum + 1721058.5:
jd0 = datenum(base.epoch) + 1721058.5; % jd0: Julian date at epoch [days]

% Mean Sun right ascension from the mean longitude of the Sun, L = 280.460 + 36000.771 T deg:
raSun = 280.460 + 36000.771*(jd0 - 2451545.0)/36525; % raSun: mean Sun right ascension at epoch [deg]

for j = 1:numel(idx)
    k = idx(j); % k: case index
    cfg = base; % cfg: parameter struct for this case
    a = base.Re + (hp(k) + ha(k))/2; % a: semi-major axis [km]
    cfg.h = a - base.Re;

    % Eccentricity from the apsis radii, e = (ra - rp)/(ra + rp):
    cfg.ecc = (ha(k) - hp(k)) / (2*a);

    % RAAN from LTAN, RAAN = RA_sun + 15 deg/h * (LTAN - 12 h):
    cfg.raan0 = mod(raSun + 15*(ltan(k) - 12), 360);

    if isnan(incFix(k))
        cfg.inc = ssoInclination(cfg);
    else
        cfg.inc = incFix(k);
    end

    % Kepler's third law:
    T = 2*pi*sqrt(a^3/base.mu); % T: orbital period [s]

    tic;
    out = propagate(cfg, nDays*86400/T); % out: state history struct from propagate.m
    ecl = eclipseTimes(out, cfg); % ecl: shadow entry and exit times and durations [s]

    res(j).name = name{k};
    res(j).label = sprintf('%s: %g x %g km, i = %.3f deg, LTAN %02d:%02d', name{k}, hp(k), ha(k), ...
                           cfg.inc, floor(ltan(k)), round(60*mod(ltan(k), 1)));
    res(j).cfg = cfg;
    res(j).T = out.T;
    res(j).tEnd = out.t(end);
    res(j).ecl = ecl;

    fprintf('%s done in %.0f s\n', res(j).label, toc);
end

end


function inc = ssoInclination(cfg)
% ssoInclination: inclination whose integrated nodal rate matches the mean motion of the Sun.
% cfg: parameter struct from config.m with the case semi-major axis, eccentricity and RAAN set
% inc: Sun synchronous inclination [deg]

a = cfg.Re + cfg.h; % a: semi-major axis [km]
p = a * (1 - cfg.ecc^2); % p: semilatus rectum [km]
n = sqrt(cfg.mu / a^3); % n: mean motion [rad/s]
T = 2*pi/n; % T: orbital period [s]
wSun = 2*pi / (365.2422*86400); % wSun: mean motion of the Sun along the ecliptic [rad/s]

% First order J2 nodal rate, dRAAN/dt = -1.5 n J2 (Re/p)^2 cos(i):
cfg.inc = acosd(-wSun / (1.5*n*cfg.J2*(cfg.Re/p)^2));

cfg.use.eclipse = false;
out = propagate(cfg, round(5*86400/T)); % out: five day calibration run from propagate.m

% RAAN from the angular momentum vector, h = |h| * [sin(i)*sin(RAAN); -sin(i)*cos(RAAN); cos(i)]:
raan = unwrap(atan2(out.hvec(1,:), -out.hvec(2,:))); % raan: 1 x N unwrapped RAAN [rad]

w = round(T/cfg.ode.dtOut); % w: one orbital period [samples]
i1 = 1:w; % i1: samples of the first orbit
i2 = numel(out.t)-w+1 : numel(out.t); % i2: samples of the last orbit

% Integrated nodal rate from orbit averaged RAAN at the start and end of the run:
rate = (mean(raan(i2)) - mean(raan(i1))) / (mean(out.t(i2)) - mean(out.t(i1))); % rate: nodal rate [rad/s]

% Nodal rate proportional to cos(i), rescaled to the Sun rate:
inc = acosd(cosd(cfg.inc) * wSun / rate);

end
