function srpCheck(nDay)
% srpCheck: validation of the solar radiation pressure in eom.m against an independent SI calculation in full sunlight, penumbra and umbra, its size against drag at the three ECSS levels, and its effect on satellites 1.1 and 2.1 over nDay days. Declared before the first run: test states at sim start, limit 1e-6 on the relative difference as for the drag in dragCheck.m, umbra exactly zero; size and effect are comparisons only.
% nDay: length of the effect runs, 30 when omitted [days]

if nargin < 1
    nDay = 30;
end

here = fileparts(mfilename('fullpath')); % here: folder holding this file
root = fullfile(here, '..'); % root: repository root folder
addpath(root);
odinPath;

outDir = fullfile(here, 'results'); % outDir: folder for logs and results
if ~isfolder(outDir)
    mkdir(outDir);
end
logFile = fullfile(outDir, sprintf('srpCheck_%dd.txt', nDay)); % logFile: text log of this run
if isfile(logFile)
    delete(logFile);
end
diary(logFile);

base = config(); % base: parameter struct from config.m
lim = 1e-6; % lim: limit on the relative difference from the SI calculation
ok = true(1, 3); % ok: pass flag of checks 1 to 3

c0 = base; % c0: gravity only, drag and solar radiation pressure off
c0.use.drag = false;
c0.use.srp = false;
c1 = c0; % c1: the same with solar radiation pressure on
c1.use.srp = true;

[rs, s0] = sunVector(0, base); % rs: Sun position at sim start [km]; s0: unit vector from the Earth centre to the Sun
rO = base.Re + base.h; % rO: orbit radius at sim start [km]
pN = cross(s0, [0; 0; 1]); % pN: direction normal to the Sun line
pN = pN/norm(pN);
vC = sqrt(base.mu/rO); % vC: circular speed [km/s]

% Points at rO on the circle through the anti-Sun point, theta from the anti-Sun direction:
rAt = @(th) rO*(-cosd(th)*s0 + sind(th)*pN); % rAt: position at angle th [km]

% Penumbra point where half the Sun disk shows, by bisection on the independent visible share:
th = [0 89]; % th: bracket, umbra at 0 and full sunlight at 89 [deg]
for k = 1:60
    tm = mean(th); % tm: midpoint [deg]
    if visibleShare(rAt(tm), rs, base) < 0.5
        th(1) = tm;
    else
        th(2) = tm;
    end
end
thP = mean(th); % thP: penumbra test angle [deg]

pts = {'full sunlight, under the Sun', rO*s0
       sprintf('penumbra, %.4f deg from the anti-Sun point', thP), rAt(thP)
       'umbra, the anti-Sun point', -rO*s0}; % pts: label and position of each test state [km]

fprintf('solar radiation pressure, spherical model: Cr %g, %g m^2 facing the Sun, %g kg, S0 %g W/m^2, c %.0f m/s\n', base.srp.Cr, base.srp.area, ...
        base.drag.mass, base.S0, base.c);
fprintf('1 to 3. eom.m with solar radiation pressure on less off, against an SI calculation with the visible share from the circle overlap, at %s:\n', ...
        datestr(datenum(base.epoch), 'dd mmm yyyy HH:MM'));
aMag = zeros(1, 3); % aMag: size of the model acceleration at each test state [m/s^2]
for j = 1:3
    r = pts{j, 2}; % r: test position [km]
    v = vC*cross([0; 0; 1], r/norm(r)); % v: a circular velocity at r, not used by the pressure [km/s]
    if norm(v) == 0
        v = [vC; 0; 0];
    end
    dA = eom(0, [r; v], c1) - eom(0, [r; v], c0); % dA: solar radiation pressure part of the state derivative [km/s; km/s^2]
    aM = 1e3*dA(4:6); % aM: model acceleration [m/s^2]
    aMag(j) = norm(aM);

    K = visibleShare(r, rs, base); % K: visible share of the Sun disk from the circle overlap
    d = rs - r; % d: spacecraft to Sun vector [km]

    % Solar radiation pressure in SI units, pressure S0/c at 1 AU scaled by the inverse square of the distance:
    aSI = -base.srp.Cr*base.srp.area/base.drag.mass*K*(base.S0/base.c)/(norm(d)/base.AU)^2*d/norm(d); % aSI: [m/s^2]

    if K == 0
        ok(j) = all(aM == 0) && all(dA(1:3) == 0);
        fprintf('   %-44s K = %.4f, model |a| %.3e m/s^2, exactly zero required: %s\n', pts{j, 1}, K, norm(aM), verdict(ok(j)));
    else
        eA = norm(aM - aSI)/norm(aSI); % eA: relative difference of the two accelerations
        ang = atan2d(norm(cross(aM, -d)), dot(aM, -d)); % ang: angle of the model acceleration from the anti-Sun line [deg]
        ok(j) = eA < lim && all(dA(1:3) == 0);
        fprintf('   %-44s K = %.4f, model |a| %.4e m/s^2, SI |a| %.4e m/s^2, relative difference %.1e, %.1e deg off the anti-Sun line, limit %.0e: %s\n', ...
                pts{j, 1}, K, norm(aM), norm(aSI), eA, ang, lim, verdict(ok(j)));
    end
end

% Size against drag around one orbit of satellite 1.1 at sim start, comparison only:
T = 2*pi*sqrt(rO^3/base.mu); % T: orbital period [s]
cP = c0; % cP: one orbit without perturbations, for the positions
cP.use.eclipse = true;
o = propagate(cP, 1); % o: one orbit of satellite 1.1 on the cfg.ode.dtOut grid
K1 = shadow(o.r, o.rsun, base); % K1: 1 x N visible share of the Sun disk around the orbit
aS = base.srp.Cr*base.srp.area*base.S0/(base.drag.mass*base.c)*K1./(vecnorm(o.rsun - o.r)/base.AU).^2; % aS: 1 x N solar radiation pressure [m/s^2]
vRel = 1e3*(o.v - base.wE*[-o.r(2, :); o.r(1, :); zeros(size(o.t))]); % vRel: 3 x N velocity relative to the co-rotating atmosphere [m/s]
F107 = [base.decay.F107Low, base.decay.F107Moderate, base.decay.F107High]; % F107: 10.7 cm flux of the ECSS low, moderate and high cases [sfu]
Ap = [base.decay.ApLow, base.decay.ApModerate, base.decay.ApHigh]; % Ap: daily geomagnetic index of the same cases
lvl = {'low', 'moderate', 'high'}; % lvl: ECSS activity level of each case
fprintf('\nsize around one orbit of satellite 1.1 from sim start (%.1f min, %.0f%% of it in shadow), orbit means, comparison only:\n', T/60, 100*mean(K1 < 1));
fprintf('   solar radiation pressure %.3e m/s^2\n', mean(aS));
aDm = zeros(1, 3); % aDm: orbit-mean drag in each case [m/s^2]
for i = 1:3
    cD = base; % cD: NRLMSISE-00 at case i
    cD.drag.model = 'NRLMSISE-00';
    cD.drag.F107 = F107(i);
    cD.drag.Ap = Ap(i);
    rho = dragDensity(o.t, o.r, cD); % rho: 1 x N density [kg/m^3]

    % Drag, |aD| = rho |vRel|^2 Cd A/(2 m):
    aDm(i) = mean(0.5*rho.*vecnorm(vRel).^2*base.drag.Cd*base.drag.area/base.drag.mass);
    fprintf('   drag at ECSS %-8s %.3e m/s^2, solar radiation pressure over drag %.2f\n', lvl{i}, aDm(i), mean(aS)/aDm(i));
end

% Effect over nDay days on satellites 1.1 and 2.1, drag off, comparison only:
fprintf('\neffect over %d days, drag off, solar radiation pressure on less off, comparison only:\n', nDay);
fprintf('   satellite   LTAN    largest position change [km]   orbit-mean altitude change [m]   largest eccentricity change   inclination change [deg]   node change [deg]   run on/off [s]\n');
eff = zeros(base.const.nPlanes, 5); % eff: position, altitude, eccentricity, inclination and node changes of each plane
for k = 1:base.const.nPlanes
    a0 = c0; % a0: satellite k.1, solar radiation pressure off
    a0.raan0 = mod(base.raan0 + (k - 1)*base.const.dRaan, 360);
    a0.use.eclipse = false;
    a1 = a0; % a1: the same with solar radiation pressure on
    a1.use.srp = true;
    nOrb = nDay*86400/T; % nOrb: run length in orbital periods

    tic;
    q0 = propagate(a0, nOrb); % q0: without solar radiation pressure
    t0 = toc; % t0: wall time without [s]
    tic;
    q1 = propagate(a1, nOrb); % q1: with solar radiation pressure
    t1 = toc; % t1: wall time with [s]

    w = round(T/base.ode.dtOut); % w: one orbital period [samples]
    alt0 = movmean(vecnorm(q0.r) - base.Re, w); % alt0: orbit-mean altitude without [km]
    alt1 = movmean(vecnorm(q1.r) - base.Re, w); % alt1: orbit-mean altitude with [km]
    kE = numel(q0.t) - w; % kE: sample one orbit before the end
    [e0, i0, W0] = elements(q0, base.mu); % e0, i0, W0: eccentricity, inclination [deg] and node [deg] without
    [e1, i1, W1] = elements(q1, base.mu); % e1, i1, W1: the same with
    iM0 = movmean(i0, w); % iM0: orbit-mean inclination without [deg]
    iM1 = movmean(i1, w); % iM1: orbit-mean inclination with [deg]
    eff(k, :) = [max(vecnorm(q1.r - q0.r)), 1e3*((alt1(kE) - alt1(w)) - (alt0(kE) - alt0(w))), max(abs(e1 - e0)), ...
                 iM1(kE) - iM0(kE), mod(W1(kE) - W0(kE) + 180, 360) - 180];
    ltan = mod(12 + (a0.raan0 - (280.460 + 36000.771*(datenum(base.epoch) + 1721058.5 - 2451545.0)/36525))/15, 24); % ltan: LTAN at sim start, mean Sun as in tradeLighting.m [h]
    mL = mod(round(60*ltan), 1440); % mL: LTAN to the nearest minute [min]
    fprintf('   %d.1         %02d:%02d   %28.3f   %30.2f   %27.2e   %24.2e   %17.2e   %.0f / %.0f\n', k, floor(mL/60), mod(mL, 60), eff(k, :), t1, t0);
end
fprintf('\nchecks passed: %d of 3\n', nnz(ok));

diary off;
save(fullfile(outDir, sprintf('srpCheck_%dd.mat', nDay)), 'ok', 'aMag', 'thP', 'aS', 'aDm', 'eff');

end


function K = visibleShare(r, rsun, cfg)
% visibleShare: share of the Sun disk not covered by the Earth disk, from the overlap area of two circles, written apart from shadow.m.
% r: satellite position in ECI [km]
% rsun: Sun position in ECI [km]
% cfg: parameter struct from config.m
% K: visible share of the Sun disk, 0 to 1

d = rsun - r; % d: satellite to Sun vector [km]
aS = asin(cfg.Rsun/norm(d)); % aS: apparent radius of the Sun [rad]
aE = asin(cfg.Re/norm(r)); % aE: apparent radius of the Earth [rad]
c = acos(max(min(dot(-r, d)/(norm(r)*norm(d)), 1), -1)); % c: angle between the Earth centre and the Sun centre [rad]

if c >= aS + aE
    K = 1;
elseif c <= aE - aS
    K = 0;
elseif c <= aS - aE
    K = 1 - (aE/aS)^2;
else
    % Lens area of two overlapping circles of radii aS and aE whose centres are c apart:
    L = aS^2*acos((c^2 + aS^2 - aE^2)/(2*c*aS)) + aE^2*acos((c^2 + aE^2 - aS^2)/(2*c*aE)) ...
        - 0.5*sqrt((-c + aS + aE)*(c + aS - aE)*(c - aS + aE)*(c + aS + aE)); % L: overlap area [rad^2]
    K = 1 - L/(pi*aS^2);
end

end


function [e, inc, raan] = elements(out, mu)
% elements: osculating eccentricity, inclination and node along a state history.
% out: state history struct from propagate.m
% mu: gravitational parameter [km^3/s^2]
% e: 1 x N eccentricity
% inc: 1 x N inclination [deg]
% raan: 1 x N right ascension of the ascending node, unwrapped [deg]

h = out.hvec; % h: 3 x N specific angular momentum [km^2/s]

% Eccentricity vector, e = v x h/mu - r/|r|:
ev = cross(out.v, h)/mu - out.r./vecnorm(out.r); % ev: 3 x N eccentricity vector
e = vecnorm(ev);
inc = acosd(h(3, :)./vecnorm(h));
raan = rad2deg(unwrap(atan2(h(1, :), -h(2, :))));

end


function s = verdict(pass)
% verdict: PASS or CHECK for a logical test result.
% pass: true when the check passed
% s: text for the log

if pass
    s = 'PASS';
else
    s = 'CHECK';
end

end
