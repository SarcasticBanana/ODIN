% UNCERTAINTY: Monte Carlo of the orbit over the mission from the Inputs sheet, insertion, mass, drag coefficient and density drawn, at ECSS low, moderate and high solar activity: altitude after the mission, LTAN, longest shadow and lowest charge, along-track error of a prediction and drift from the planned slot, and time in view of the ground stations a day, each shortcut checked against ode45 or MISSION.m with pass marks set before the first run. Figures in results/uncertainty at 14 pt, every number in results/uncertainty.txt and results/uncertainty.json, the runs in results/uncertainty.mat. Needs results/mission.mat from MISSION.m.

clear; close all; clc;

here = fileparts(mfilename('fullpath')); % here: repository root folder
addpath(here);
odinPath;

rerun = true; % rerun: false redraws the figures and rewrites the numbers from results/uncertainty.mat without the runs
hLine = 550; % hLine: altitude the five-year decay is judged against, the lowest the battery is sized for (cfg.batt.h), as in MISSION.m [km]
comMin = 190; % comMin: communications time a day the science model takes from Comms, DR2 slide 41 [min]
dAhead = [1 3 7]; % dAhead: horizons of the along-track error of a prediction [days]
pct = [0.135 50 99.865]; % pct: percentiles reported, the median between the 3-sigma ends of a normal distribution [%]
p3 = 99.73; % p3: percentile of a magnitude reported as its 3-sigma value [%]
fs = 14; % fs: font size of every figure, the 14 pt floor of the DR2 assignment [pt]

hT = {540:20:660, 400:20:660, 180:20:660}; % hT: heights of the drag table at ECSS low, moderate and high, semi-major axis less cfg.Re [km]
dT = 20; % dT: date step of the drag table [days]
nUT = 4; % nUT: start times averaged at each date of the drag table, spread evenly over the day
lT = {-3:3, -3:5, -3:12}; % lT: LTAN offsets of the drag table at ECSS low, moderate and high [h]
bP = 0:0.05:90; % bP: Sun angles above the orbit plane in the shadow table [deg]
hP = 180:5:660; % hP: heights of the shadow table [km]
hC = [200 250 300 350 400 450 500 550 575 600 615 625 635 650]; % hC: altitudes of the contact runs, osculating at sim start [km]
nDayC = 10; % nDayC: length of each contact run [days]
chkInc = 0.1; % chkInc: inclination step of the node rate check [deg]
chkH = 10; % chkH: altitude step of the node rate check [km]
nDayW = 20; % nDayW: length of the node rate runs [days]
nPath = 50; % nPath: ode45 runs against the along-track error of the GNSS fix
nDayP = 3; % nDayP: length of the along-track runs [days]
dragX = 1.5; % dragX: drag coefficient of the second along-track drag run over the sheet's

% Pass marks of the checks, set before the first run:
tolDecayKm = 0.2; % tolDecayKm: table decay of the sheet spacecraft against dragDecay.m after the mission [km]
tolDecayRel = 0.01; % tolDecayRel: or this share of the drop, whichever is larger
tolNode = 0.02; % tolNode: node rate change against ode45, share
tolShadowMin = 0.1; % tolShadowMin: longest shadow of each year against MISSION.m [min]
tolShadowDay = 3; % tolShadowDay: and its day [days]
tolPath = 0.10; % tolPath: along-track against ode45, the GNSS rms difference over the rms spread and the drag part, share
tolContact = 0.01; % tolContact: contact time change at the 3-sigma inclinations, share
tolConv = 0.02; % tolConv: change of the 3-sigma percentiles from half the draws to all, share of the 3-sigma spread

outDir = fullfile(here, 'results'); % outDir: folder for the numbers and the saved runs
figDir = fullfile(outDir, 'uncertainty'); % figDir: folder for the slide figures
if ~isfolder(figDir)
    mkdir(figDir);
end
matFile = fullfile(outDir, 'uncertainty.mat'); % matFile: saved runs

col = [0.00 0.45 0.70
       0.84 0.37 0.00
       0.00 0.62 0.45
       0.80 0.47 0.65
       0.90 0.62 0.00
       0.34 0.71 0.91]; % col: 6 x 3 plane colours, as in MISSION.m
caseCol = col([3 5 2], :); % caseCol: 3 x 3 colour of the ECSS low, moderate and high cases, as in MISSION.m
caseName = {'low', 'moderate', 'high'}; % caseName: ECSS activity level of each case

base = config(); % base: parameter struct of ODIN 1.1 from the Inputs sheet
mc = base.mc; % mc: Monte Carlo inputs from the sheet
st = stations(); % st: ground stations in use
M = load(fullfile(outDir, 'mission.mat'), 'sv', 'dec', 'b', 'cB'); % M: MISSION.m runs, shadow surveys, decays, battery and the parameter struct they ran with

% MISSION.m must have run with this sheet's orbit, spacecraft and activity levels:
same = [M.cB.h M.cB.inc M.cB.raan0 M.cB.drag.mass M.cB.drag.Cd M.cB.drag.area M.cB.decay.F107Low M.cB.decay.F107Moderate M.cB.decay.F107High ...
        M.cB.decay.ApLow M.cB.decay.ApModerate M.cB.decay.ApHigh M.cB.power.years M.cB.decay.dt];
now0 = [base.h base.inc base.raan0 base.drag.mass base.drag.Cd base.drag.area base.decay.F107Low base.decay.F107Moderate base.decay.F107High ...
        base.decay.ApLow base.decay.ApModerate base.decay.ApHigh base.power.years base.decay.dt];
if any(abs(same - now0) > 1e-9)
    error('UNCERTAINTY:mission', 'results/mission.mat was run with other inputs than the sheet, run MISSION.m first');
end

yrs = base.power.years; % yrs: mission length [yr]
nDay = round(365.25*yrs); % nDay: days in the mission
nPl = base.const.nPlanes; % nPl: number of orbit planes
dt = 86400*base.decay.dt; % dt: decay step [s]
nStep = ceil(365.25*yrs*86400/dt); % nStep: decay steps over the mission, as in dragDecay.m
kDay = round(86400/dt); % kDay: decay steps a day
kBand = 10*kDay; % kBand: decay steps between the stored percentiles, 10 days
kYear = round(365.25*86400*(1:yrs)/dt); % kYear: decay step ending each mission year
dn0 = datenum(base.epoch); % dn0: date number at sim start [days]
F107 = [base.decay.F107Low, base.decay.F107Moderate, base.decay.F107High]; % F107: 10.7 cm solar flux of the ECSS low, moderate and high cases [sfu]
Ap = [base.decay.ApLow, base.decay.ApModerate, base.decay.ApHigh]; % Ap: daily geomagnetic index of the same cases
wSun = 36000.771/36525; % wSun: mean Sun rate along the equator, as in the LTAN of MISSION.m [deg/day]
h0 = M.sv{1}.h0; % h0: orbit-mean altitude of plane 1 at sim start, the start of the MISSION.m decays [km]
iS = size(M.dec, 2) - 1; % iS: column of M.dec holding the sheet spacecraft
F0 = base.drag.Cd*base.drag.area/base.drag.mass; % F0: Cd A/m of the sheet spacecraft [m^2/kg]
raanP = mod(base.raan0 + (0:nPl-1)*base.const.dRaan, 360); % raanP: node of each plane at sim start, as in constellation.m [deg]
hC = unique([hC, base.h]); % hC: contact altitudes with the sheet's added, the reference of the inclination check [km]
eFullE = base.batt.energy*base.batt.socTop*M.b.qEOL; % eFullE: full charge at the end of life, as in MISSION.m [Wh]

% Secular J2 node rate of a circular orbit, -(3/2) J2 (Re/a)^2 n cos i, as in dragDecay.m:
wdot = @(a, inc) -1.5*base.J2*(base.Re./a).^2.*sqrt(base.mu./a.^3).*cosd(inc)*180/pi*86400; % wdot: node rate [deg/day] from a [km] and i [deg]

% Secular rate of the argument of latitude, n (1 + (3/2) J2 (Re/a)^2 (4 cos^2 i - 1)), as in dragDecay.m:
udot = @(a, inc) sqrt(base.mu./a.^3).*(1 + 1.5*base.J2*(base.Re./a).^2.*(4*cosd(inc).^2 - 1)); % udot: [rad/s] from a [km] and i [deg]

% Mean semi-major axis whose secular node rate is that of the ode45 orbit, so the draws are taken as differences from it:
aRef = fzero(@(a) wdot(a, base.inc) - M.sv{1}.dW, base.Re + base.h); % aRef: [km]
nRef = udot(aRef, base.inc); % nRef: rate of the argument of latitude of the reference orbit [rad/s]

% ---------------------------------------------------------------- draws
rng(mc.seed, 'twister');
N = mc.n; % N: number of draws
N1 = N + 1; % N1: draws plus the sheet spacecraft, last
dHp = mc.perigee3/3*randn(N, 1); % dHp: perigee insertion error [km]
dHa = mc.apogee3/3*randn(N, 1); % dHa: apogee insertion error [km]
dInc = mc.inc3/3*randn(N, 1); % dInc: inclination insertion error [deg]
dRaan = mc.raan3/3*randn(N, 1); % dRaan: RAAN insertion error, instantaneous launch window [deg]
mass = mc.massMin + (base.drag.mass - mc.massMin)*rand(N, 1); % mass: spacecraft mass [kg]
Cd = mc.CdMin + (base.drag.Cd - mc.CdMin)*rand(N, 1); % Cd: drag coefficient
fRho = exp(mc.rhoSigma*randn(N, 1)); % fRho: density over NRLMSISE-00, lognormal
gn = [mc.gnssPos*randn(N, 3), mc.gnssVel*randn(N, 3)]; % gn: GNSS fix error, radial, along-track and cross-track position [m] and velocity [m/s]
dRaanW = mc.raan3Window/3*randn(N, 1); % dRaanW: RAAN insertion error with a 1-hour launch window [deg]

% Circular orbit at the mean of perigee and apogee:
dh = (dHp + dHa)/2; % dh: altitude error [km]

F = [fRho.*Cd*base.drag.area./mass; F0]; % F: N1 x 1 Cd A/m times the density factor [m^2/kg]
hS0 = h0 + [dh; 0]; % hS0: N1 x 1 orbit-mean altitude at sim start [km]
inc = base.inc + [dInc; 0]; % inc: N1 x 1 inclination [deg]
L0 = [dRaan; 0]/15; % L0: N1 x 1 LTAN offset at sim start, 15 deg/h [h]
dWW = [dRaanW - dRaan; 0]; % dWW: N1 x 1 node of the 1-hour window case less that of the first case [deg]

% LTAN rate of a draw, its secular node rate less that of the reference orbit, 15 deg/h:
rateL = @(h, in) (wdot(aRef + h - h0, in) - wdot(aRef, base.inc))/15; % rateL: [h/day] from the orbit-mean altitude [km] and inclination [deg]

% LTAN rate of the node of dragDecay.m, its secular rate at a = Re + h less the mean Sun:
rateD = @(h, in) (wdot(base.Re + h, in) - wSun)/15; % rateD: [h/day]

tb = [0:kBand:nStep-1, nStep]*dt/86400; % tb: times of the stored percentiles [days]

if rerun
    tic;

    % ------------------------------------------------------------ tables
    tT = 0:dT:dT*ceil(365.25*yrs/dT + 1); % tT: dates of the drag table, past the last midpoint [days]
    tabs = cell(1, 3); % tabs: drag table of each case
    for i = 1:3
        c = base; % c: parameter struct of case i
        c.drag.F107 = F107(i);
        c.drag.Ap = Ap(i);
        tabs{i} = rateTable(c, hT{i}, lT{i}, tT, nUT, wSun);
        fprintf('drag table at ECSS %s, %d orbits, %.0f s so far\n', caseName{i}, numel(hT{i})*numel(lT{i})*numel(tT)*nUT, toc);
    end
    phi = shadowTable(base, bP, hP); % phi: half-width of the shadow arc against Sun angle and height

    % ------------------------------------------------------------ the draws over the mission
    H5 = NaN(N1, 3); % H5: altitude at the end of the mission, or where it fell below the floor [km]
    LE = NaN(N1, 3); % LE: LTAN offset at the end [h]
    T6 = NaN(N1, 3); % T6: time below hLine, NaN if never [s]
    TF = NaN(N1, 3); % TF: time below cfg.decay.hFloor, NaN if never [s]
    HY = NaN(N1, yrs, 3); % HY: altitude at the end of each mission year [km]
    LY = NaN(N1, yrs, 3); % LY: LTAN offset at the end of each mission year [h]
    DM = zeros(N1, 3); % DM: longest shadow of either plane over the mission [s]
    DW = zeros(N1, 3); % DW: the same with the 1-hour window RAAN error [s]
    DD = NaN(N1, 3); % DD: day of the longest shadow, 0 on the first day
    DB = NaN(N1, 3); % DB: Sun angle above the orbit plane there [deg]
    DP = NaN(N1, 3); % DP: plane of the longest shadow
    DH = NaN(N1, 3); % DH: altitude there [km]
    BH = NaN(3, numel(tb), 3); % BH: percentiles pct of the altitude of the draws at tb, each case [km]
    BL = NaN(3, numel(tb), 3); % BL: percentiles pct of the LTAN offset at tb [h]
    nOut = zeros(1, 3); % nOut: draw steps outside the drag table, each case
    for i = 1:3
        hS = hS0; % hS: N1 x 1 orbit-mean altitude [km]
        L = L0; % L: N1 x 1 LTAN offset [h]
        tF = NaN(N1, 1); % tF: N1 x 1 time below the floor [s]
        t6 = NaN(N1, 1); % t6: N1 x 1 time below hLine [s]
        t6(hS0 < hLine) = 0;
        jb = 1; % jb: next column of the stored percentiles
        for k = 0:nStep-1
            t = k*dt; % t: time at the start of the step [s]

            if mod(k, kBand) == 0
                BH(:, jb, i) = pctl(hS(1:N), pct);
                BL(:, jb, i) = pctl(L(1:N), pct);
                jb = jb + 1;
            end

            % Longest shadow of the day, at 00:00 UTC, each plane:
            if mod(k, kDay) == 0 && k/kDay < nDay
                [~, s] = sunVector(t, base); % s: unit vector to the Sun
                q = find(isnan(tF)); % q: draws above the floor
                for p = 1:nPl
                    W = raanP(p) + M.sv{p}.dW*t/86400 + 15*L(q); % W: node of plane p [deg]
                    uD = 2*pi/M.sv{p}.TN*udot(aRef + hS(q) - h0, inc(q))/udot(aRef, base.inc); % uD: rate of the argument of latitude, ode45 nodal period at the start [rad/s]
                    [d1, b1] = shadowLength(s, hS(q), inc(q), W, uD, phi); % d1: shadow length [s], b1: Sun angle [deg]
                    up = d1 > DM(q, i); % up: draws with a new longest shadow
                    DM(q(up), i) = d1(up);
                    DD(q(up), i) = k/kDay;
                    DB(q(up), i) = b1(up);
                    DP(q(up), i) = p;
                    DH(q(up), i) = hS(q(up));
                    d2 = shadowLength(s, hS(q), inc(q), W + dWW(q), uD, phi); % d2: the 1-hour window case [s]
                    DW(q, i) = max(DW(q, i), d2);
                end
            end

            hB = hS; % hB: altitude at the start of the step [km]
            [hS, L, tF, o] = decayStep(hS, L, tF, t, F, inc, tabs{i}, rateL, base);
            nOut(i) = nOut(i) + o;

            x = isnan(t6) & hS < hLine; % x: draws crossing hLine in this step
            t6(x) = t + (hLine - hB(x))./(hS(x) - hB(x))*dt;

            y = find(kYear == k + 1); % y: mission year ending with this step
            if ~isempty(y)
                HY(:, y, i) = hS;
                LY(:, y, i) = L;
            end
        end
        BH(:, end, i) = pctl(hS(1:N), pct);
        BL(:, end, i) = pctl(L(1:N), pct);
        H5(:, i) = hS;
        LE(:, i) = L;
        T6(:, i) = t6;
        TF(:, i) = tF;
        fprintf('ECSS %s: %d draws over %g years, %.0f s so far\n', caseName{i}, N, yrs, toc);
    end

    % ------------------------------------------------------------ check 1: table decay of the sheet spacecraft against dragDecay.m
    hK = NaN(3, nStep + 1); % hK: table decay of the sheet spacecraft with the node of dragDecay.m [km]
    for i = 1:3
        hs = h0; % hs: altitude [km]
        l = 0; % l: LTAN offset [h]
        tf = NaN; % tf: time below the floor [s]
        hK(i, 1) = hs;
        for k = 0:nStep-1
            [hs, l, tf] = decayStep(hs, l, tf, k*dt, F0, base.inc, tabs{i}, rateD, base);
            hK(i, k + 2) = hs;
        end
    end

    % ------------------------------------------------------------ check 2: node rate change against ode45
    cW = base; % cW: parameter struct of the node rate runs, drag off
    cW.use.drag = false;
    wRun = zeros(1, 3); % wRun: node rate of each run, nominal, inclination step, altitude step [deg/day]
    for j = 1:3
        c = cW;
        if j == 2
            c.inc = base.inc + chkInc;
        elseif j == 3
            c.h = base.h + chkH;
        end
        [r0, v0] = oe2rv(c.Re + c.h, c.ecc, c.inc, c.raan0, c.argp0, c.nu0, c.mu);
        % r0: initial position [km], v0: initial velocity [km/s]
        [tw, Xw] = odeRun(c, [r0; v0], nDayW); % tw: output times [s], Xw: states [km, km/s]
        hw = cross(Xw(:, 1:3), Xw(:, 4:6), 2); % hw: angular momentum [km^2/s]
        pw = polyfit(tw/86400, unwrap(atan2(hw(:, 1), -hw(:, 2)))*180/pi, 1); % pw: node against days, slope [deg/day]
        wRun(j) = pw(1);
    end
    fprintf('node rate runs done, %.0f s so far\n', toc);

    % ------------------------------------------------------------ check 3: nominal shadow, drag off, against MISSION.m
    DS = zeros(nPl, nDay); % DS: semi-analytic longest shadow of each day of the drag-free sheet orbit [s]
    for d = 0:nDay-1
        t = d*86400;
        [~, s] = sunVector(t, base);
        for p = 1:nPl
            DS(p, d + 1) = shadowLength(s, h0, base.inc, raanP(p) + M.sv{p}.dW*t/86400, 2*pi/M.sv{p}.TN, phi);
        end
    end

    % ------------------------------------------------------------ check 4: along-track of the GNSS fix against ode45
    cP = base; % cP: parameter struct of the along-track runs, drag off
    cP.use.drag = false;
    [r0, v0] = oe2rv(base.Re + base.h, base.ecc, base.inc, base.raan0, base.argp0, base.nu0, base.mu);
    eR = r0/norm(r0); % eR: radial unit vector at sim start
    eC = cross(r0, v0)/norm(cross(r0, v0)); % eC: cross-track unit vector
    eA = cross(eC, eR); % eA: along-track unit vector
    [tp, X0] = odeRun(cP, [r0; v0], nDayP); % tp: output times [s], X0: nominal states
    kP = round((1:nDayP)*86400/base.ode.dtOut) + 1; % kP: output samples at whole days
    yOde = zeros(nPath, nDayP); % yOde: along-track of each run ahead of the nominal at whole days [m]
    for j = 1:nPath
        dr = [eR eA eC]*gn(j, 1:3).'/1e3; % dr: position error [km]
        dv = [eR eA eC]*gn(j, 4:6).'/1e3; % dv: velocity error [km/s]
        [~, Xj] = odeRun(cP, [r0 + dr; v0 + dv], nDayP);
        yOde(j, :) = alongTrack(X0(kP, :), Xj(kP, :));
    end
    fprintf('along-track runs done, %.0f s so far\n', toc);

    % ------------------------------------------------------------ check 5: along-track of the drag part against ode45, ECSS moderate
    cG = base; % cG: parameter struct of the drag runs
    cG.use.drag = true;
    cG.drag.F107 = F107(2);
    cG.drag.Ap = Ap(2);
    [~, XG1] = odeRun(cG, [r0; v0], nDayP); % XG1: states with the sheet's drag
    cG.drag.Cd = dragX*base.drag.Cd;
    [~, XG2] = odeRun(cG, [r0; v0], nDayP); % XG2: states with dragX times the drag
    yDragOde = alongTrack(XG1(kP, :), XG2(kP, :)); % yDragOde: along-track of the stronger drag run ahead [m]
    fprintf('drag runs done, %.0f s so far\n', toc);

    % ------------------------------------------------------------ contacts
    cMin = zeros(1, numel(hC)); % cMin: time in view of any station a day, each altitude [min]
    hMean = zeros(1, numel(hC)); % hMean: orbit-mean altitude of each run [km]
    cDay = zeros(numel(hC), nDayC); % cDay: time in view each day [min]
    for j = 1:numel(hC)
        c = base;
        c.use.drag = false;
        c.h = hC(j);
        [cMin(j), hMean(j), cDay(j, :)] = contactRun(c, nDayC, st);
    end
    cInc = zeros(1, 2); % cInc: time in view a day at the sheet altitude with the 3-sigma inclinations [min]
    for j = 1:2
        c = base;
        c.use.drag = false;
        c.inc = base.inc + (2*j - 3)*mc.inc3;
        cInc(j) = contactRun(c, nDayC, st);
    end
    fprintf('contact runs done, %.0f s in all\n', toc);

    save(matFile, 'tabs', 'phi', 'H5', 'LE', 'T6', 'TF', 'HY', 'LY', 'DM', 'DW', 'DD', 'DB', 'DP', 'DH', 'BH', 'BL', 'nOut', 'hK', 'wRun', 'DS', ...
         'yOde', 'yDragOde', 'cMin', 'hMean', 'cDay', 'cInc', '-v7');
else
    load(matFile);
end

% ---------------------------------------------------------------- numbers
R = struct();
R.version = '1.1.13';
R.run = datestr(now, 'yyyy-mm-dd HH:MM');
R.matlab = version;
txtFile = fullfile(outDir, 'uncertainty.txt'); % txtFile: printed summary
if isfile(txtFile)
    delete(txtFile);
end
diary(txtFile);
fprintf('UNCERTAINTY, ODIN-Model %s, run %s, MATLAB %s\n\n', R.version, R.run, R.matlab);

R.inputs.n = N;
R.inputs.seed = mc.seed;
R.inputs.alt3 = 3*std(dh); % 3-sigma altitude error of the draws [km]
R.inputs.alt3Exact = sqrt((mc.perigee3/3)^2 + (mc.apogee3/3)^2)/2*3; % the same from the sheet [km]
R.inputs.inc3 = mc.inc3;
R.inputs.raan3 = mc.raan3;
R.inputs.raan3Window = mc.raan3Window;
R.inputs.mass = [mc.massMin base.drag.mass];
R.inputs.Cd = [mc.CdMin base.drag.Cd];
R.inputs.area = base.drag.area;
R.inputs.rhoSigma = mc.rhoSigma;
R.inputs.gnss = [mc.gnssPos mc.gnssVel];
R.inputs.B = [base.drag.mass/(base.drag.Cd*base.drag.area), pctl(1./F(1:N), pct)]; % m/(Cd A) over the density factor, sheet then percentiles [kg/m^2]
R.inputs.F = [F0, pctl(F(1:N), pct)]/F0; % drag of the draws over the sheet's
R.inputs.perigee3 = mc.perigee3;
R.inputs.apogee3 = mc.apogee3;
R.inputs.sepVel = mc.sepVel;
R.limits = struct('hLine', hLine, 'tEcl', base.batt.tEcl, 'socSafe', base.modes.socSafe, 'comMin', comMin, 'years', yrs, 'dU', base.const.dU, ...
                  'hFloor', base.decay.hFloor); % limits the results are judged against
fprintf('%d draws, seed %d; altitude 3 sigma %.2f km (perigee %g km, apogee %g km, mean of the two), inclination %g deg, RAAN %g deg (%g deg with a 1-hour window)\n', ...
        N, mc.seed, R.inputs.alt3, mc.perigee3, mc.apogee3, mc.inc3, mc.raan3, mc.raan3Window);
fprintf('mass %g to %g kg, Cd %g to %g, area %g m^2, density factor lognormal sigma %g: drag over the sheet''s %.2f, %.2f, %.2f at %s percent\n', ...
        mc.massMin, base.drag.mass, mc.CdMin, base.drag.Cd, base.drag.area, mc.rhoSigma, R.inputs.F(2:4), mat2str(pct));
fprintf('GNSS fix sigma %g m and %g m/s each axis\n\n', mc.gnssPos, mc.gnssVel);

sel = 1:N; % sel: the draws, the sheet spacecraft left out
socOf = @(d) 1 - base.modes.pNominal*d/3600/M.b.etaDisEOL/eFullE; % socOf: lowest charge at the end of life through a shadow of d seconds at the Nominal load, as in MISSION.m

% Contacts at each draw's altitude, the reference runs splined, none once below the floor:
cAt = @(h) interp1(hMean, cMin, h, 'pchip', 'extrap'); % cAt: time in view a day [min] at an orbit-mean altitude [km]
cS = cAt(hS0); % cS: N1 x 1 time in view a day at sim start [min]
R.contact.hMean = hMean;
R.contact.minPerDay = cMin;
R.contact.dayMin = min(cDay, [], 2).';
R.contact.dayMax = max(cDay, [], 2).';
R.contact.start = [cS(N1), pctl(cS(sel), pct)];
R.contact.comMin = comMin;

fprintf('altitude after %g years, time below %g km, below the %g km floor, LTAN offset, longest shadow and lowest end-of-life charge, nominal then %s percent:\n', ...
        yrs, hLine, base.decay.hFloor, mat2str(pct));
for i = 1:3
    C = struct(); % C: numbers of case i
    C.level = caseName{i};
    C.F107 = F107(i);
    C.Ap = Ap(i);
    C.h5 = [H5(N1, i), pctl(H5(sel, i), pct)]; % altitude after the mission [km]
    C.below600 = 100*mean(~isnan(T6(sel, i))); % share of draws below hLine within the mission [%]
    t6y = T6(:, i)/86400/365.25; % t6y: time below hLine [yr]
    t6y(isnan(t6y)) = Inf;
    C.t600 = [t6y(N1), pctl(t6y(sel), pct)]; % time below hLine, Inf if never [yr]
    C.fell = 100*mean(~isnan(TF(sel, i))); % share of draws below the floor within the mission [%]
    C.ltan = 60*[LE(N1, i), pctl(LE(sel, i), pct)]; % LTAN offset at the end or the floor [min]
    C.hYear = zeros(yrs, 4); % altitude at the end of each mission year, nominal then percentiles [km]
    C.ltanYear = zeros(yrs, 4); % LTAN offset at the end of each mission year [min]
    for y = 1:yrs
        C.hYear(y, :) = [HY(N1, y, i), pctl(HY(sel, y, i), pct)];
        C.ltanYear(y, :) = 60*[LY(N1, y, i), pctl(LY(sel, y, i), pct)];
    end
    C.shadow = [DM(N1, i), pctl(DM(sel, i), pct)]/60; % longest shadow [min]
    C.shadowWindow = [DW(N1, i), pctl(DW(sel, i), pct)]/60; % longest shadow with the 1-hour window [min]
    C.shadowNom = struct('day', DD(N1, i), 'date', datestr(dn0 + DD(N1, i), 'dd mmm yyyy'), 'beta', DB(N1, i), 'plane', DP(N1, i), 'h', DH(N1, i));
    [~, jw] = max(DM(sel, i)); % jw: draw with the longest shadow
    C.shadowWorst = struct('day', DD(jw, i), 'date', datestr(dn0 + DD(jw, i), 'dd mmm yyyy'), 'beta', DB(jw, i), 'plane', DP(jw, i), 'h', DH(jw, i), ...
                           'dur', DM(jw, i)/60);
    qo = DM(sel, i)/60 > base.batt.tEcl; % qo: draws with a shadow longer than the battery sizing one
    C.overSizing = 100*mean(qo); % their share [%]
    C.overSizingHmax = max([DH(qo, i); NaN]); % highest altitude of such a longest shadow [km]
    C.soc = 100*socOf([DM(N1, i), fliplr(pctl(DM(sel, i), pct))]); % lowest charge at the end of life, nominal then percentiles [%]
    C.socWindow = 100*socOf([DW(N1, i), fliplr(pctl(DW(sel, i), pct))]);
    c5 = cAt(H5(:, i)); % c5: time in view a day at the end [min]
    c5(~isnan(TF(:, i))) = 0;
    C.contact5 = [c5(N1), pctl(c5(sel), pct)];
    C.tableOut = nOut(i);
    R.levels(i) = C;
    fprintf('   %-8s %6.1f km [%6.1f %6.1f %6.1f]; below %g km: %5.1f%%, time %s yr; floor %4.1f%%; LTAN %+6.1f min [%+6.1f %+6.1f %+6.1f]; shadow %5.2f min [%5.2f %5.2f %5.2f], window 99.865: %5.2f; charge %4.1f%% [%4.1f %4.1f %4.1f]\n', ...
            C.level, C.h5, hLine, C.below600, mat2str(round(100*C.t600)/100), C.fell, C.ltan, C.shadow, C.shadowWindow(4), C.soc);
    fprintf('   %-8s time in view a day after the mission %.1f min [%.1f %.1f %.1f]; longest shadow of the draws %.2f min on %s, beta %.1f deg, %.0f km, plane %d\n', ...
            '', C.contact5, C.shadowWorst.dur, C.shadowWorst.date, C.shadowWorst.beta, C.shadowWorst.h, C.shadowWorst.plane);
    fprintf('   %-8s longest shadow above the %.1f min battery sizing shadow: %.2f%% of draws, all at or below %.0f km\n', '', base.batt.tEcl, ...
            C.overSizing, C.overSizingHmax);
end
fprintf('   battery sizing shadow %.1f min, Safe below %g%% charge\n', base.batt.tEcl, 100*base.modes.socSafe);
fprintf('time in view of any station a day at sim start: %.1f min [%.1f %.1f %.1f]; the science model takes %g min a day\n\n', R.contact.start, comMin);

% Along-track of a prediction, Clohessy-Wiltshire from the GNSS fix, rates in the frame turning with the orbit:
tA = dAhead*86400; % tA: horizons [s]
x0 = gn(:, 1); % x0: radial position error [m]
y0 = gn(:, 2); % y0: along-track position error [m]
xd = gn(:, 4) + nRef*y0; % xd: radial rate [m/s]
yd = gn(:, 5) - nRef*x0; % yd: along-track rate [m/s]
cw = @(t) 6*(sin(nRef*t) - nRef*t).*x0 + y0 - 2/nRef*(1 - cos(nRef*t)).*xd + (4*sin(nRef*t) - 3*nRef*t)/nRef.*yd; % cw: along-track at t [s], N x 1 [m]

% Drag rate difference of each draw from the sheet's at sim start, along-track -(3/4) n (da/dt) t^2:
K0 = zeros(N, 3); % K0: rate per unit drag at each draw's start altitude, LTAN offset 0 [km/s per m^2/kg]
for i = 1:3
    K0(:, i) = tableRate(tabs{i}, hS0(1:N), zeros(N, 1), 0);
end
dA = 1e3*(F(1:N) - F0).*K0; % dA: N x 3 drag rate of change of semi-major axis less the sheet's [m/s]
yDrag = @(t, i) -0.75*nRef*dA(:, i)*t.^2; % yDrag: along-track at t [s] in case i, N x 1 [m]

R.path.days = dAhead;
R.path.gnss = zeros(1, numel(dAhead)); % 3-sigma along-track of the GNSS fix [km]
R.path.drag = zeros(3, numel(dAhead)); % 3-sigma along-track of the drag part, each case [km]
R.path.total = zeros(3, numel(dAhead)); % 3-sigma along-track of both [km]
for a = 1:numel(dAhead)
    yc = cw(tA(a)); % yc: along-track of the GNSS fix [m]
    R.path.gnss(a) = pctl(abs(yc), p3)/1e3;
    for i = 1:3
        yg = yDrag(tA(a), i); % yg: along-track of the drag part [m]
        R.path.drag(i, a) = pctl(abs(yg), p3)/1e3;
        R.path.total(i, a) = pctl(abs(yc + yg), p3)/1e3;
    end
end

% Drift from the planned slot, the change of the argument of latitude rate at sim start:
drift = (udot(aRef + dh, base.inc + dInc) - nRef)*86400*180/pi; % drift: N x 1 [deg/day]
R.path.drift = [pctl(drift, pct), pctl(abs(drift), p3)];
v0n = sqrt(base.mu/aRef); % v0n: circular speed of the reference orbit [km/s]

% Separation velocity along the track, da = 2 a^2 v dv/mu, drift (3/2) n da/a:
aSep = 2*aRef^2*v0n*mc.sepVel/1e3/base.mu; % aSep: semi-major axis change [km]
R.path.sepDrift = 1.5*nRef*aSep/aRef*86400*180/pi; % [deg/day]
R.path.sepA = aSep;
fprintf('along-track of a prediction, 3 sigma [km] at %s days: GNSS fix %s\n', mat2str(dAhead), mat2str(round(100*R.path.gnss)/100));
for i = 1:3
    fprintf('   %-8s drag part %s, both %s\n', caseName{i}, mat2str(round(1000*R.path.drag(i, :))/1000), mat2str(round(100*R.path.total(i, :))/100));
end
fprintf('drift from the planned slot %.2f deg/day at 3 sigma (%s at %s percent); %g m/s of separation along the track: %.2f km, %.2f deg/day\n\n', ...
        R.path.drift(4), mat2str(round(100*R.path.drift(1:3))/100), mat2str(pct), mc.sepVel, aSep, R.path.sepDrift);

% ---------------------------------------------------------------- checks
chk = struct('name', {}, 'value', {}, 'limit', {}, 'pass', {}); % chk: every check, value against its limit

for i = 1:3
    hd = M.dec{i, iS}; % hd: dragDecay.m run of the sheet spacecraft in case i
    m = min(numel(hd.h), size(hK, 2)); % m: common steps
    e = hK(i, m) - hd.h(m); % e: table decay less dragDecay.m at the end [km]
    lim = max(tolDecayKm, tolDecayRel*(hd.h(1) - hd.h(m))); % lim: pass mark [km]
    chk(end + 1) = struct('name', sprintf('decay table against dragDecay, ECSS %s, altitude at the end [km]', caseName{i}), 'value', e, 'limit', lim, ...
                          'pass', abs(e) <= lim);
end

dSec = [wdot(aRef, base.inc + chkInc), wdot(aRef + chkH, base.inc)] - wdot(aRef, base.inc); % dSec: secular node rate change [deg/day]
dOde = wRun(2:3) - wRun(1); % dOde: ode45 node rate change [deg/day]
stepName = {sprintf('+%g deg', chkInc), sprintf('+%g km', chkH)}; % stepName: step of each run
for j = 1:2
    e = dSec(j)/dOde(j) - 1;
    chk(end + 1) = struct('name', sprintf('node rate change for %s, secular over ode45 less 1', stepName{j}), 'value', e, 'limit', tolNode, ...
                          'pass', abs(e) <= tolNode);
end

for p = 1:nPl
    eY = zeros(1, ceil(yrs)); % eY: longest shadow of each year less MISSION.m [min]
    eD = zeros(1, ceil(yrs)); % eD: its day less MISSION.m's
    for y = 1:ceil(yrs)
        qd = find((1:nDay) > (y - 1)*365.25 & (1:nDay) <= y*365.25); % qd: days of year y, as in MISSION.m
        [a1, j1] = max(DS(p, qd));
        [a2, j2] = max(M.sv{p}.day(qd));
        eY(y) = (a1 - a2)/60;
        eD(y) = j1 - j2;
    end
    [~, y1] = max(abs(eY));
    chk(end + 1) = struct('name', sprintf('longest shadow of each year, plane %d, against MISSION [min]', p), 'value', eY(y1), 'limit', tolShadowMin, ...
                          'pass', all(abs(eY) <= tolShadowMin));
    [~, y2] = max(abs(eD));
    chk(end + 1) = struct('name', sprintf('day of the longest shadow of each year, plane %d, against MISSION [days]', p), 'value', eD(y2), 'limit', tolShadowDay, ...
                          'pass', all(abs(eD) <= tolShadowDay));
end

yCw = zeros(nPath, nDayP); % yCw: Clohessy-Wiltshire along-track of the first nPath draws at whole days [m]
for d = 1:nDayP
    yc = cw(d*86400);
    yCw(:, d) = yc(1:nPath);
end
e = sqrt(mean((yOde(:, end) - yCw(:, end)).^2))/sqrt(mean(yCw(:, end).^2)); % e: rms difference over rms spread at the last day
chk(end + 1) = struct('name', sprintf('along-track of the GNSS fix after %d days, rms of ode45 less Clohessy-Wiltshire over rms, %d runs', nDayP, nPath), ...
                      'value', e, 'limit', tolPath, 'pass', e <= tolPath);

KG = 0; % KG: rate per unit drag at the start altitude, mean over the drag runs [km/s per m^2/kg]
for j = 0:12
    KG = KG + tableRate(tabs{2}, h0, 0, j*nDayP*86400/12)/13;
end
ySec = -0.75*nRef*1e3*(dragX - 1)*F0*KG*(nDayP*86400)^2; % ySec: secular along-track of the stronger drag run after nDayP days [m]
e = ySec/yDragOde(end) - 1;
chk(end + 1) = struct('name', sprintf('along-track of %g times the drag after %d days, ECSS moderate, secular over ode45 less 1', dragX, nDayP), ...
                      'value', e, 'limit', tolPath, 'pass', abs(e) <= tolPath);

jC = find(hC == base.h, 1); % jC: contact run at the sheet altitude
e = max(abs(cInc/cMin(jC) - 1));
chk(end + 1) = struct('name', sprintf('time in view a day at %g deg either side of the inclination, change', mc.inc3), 'value', e, 'limit', tolContact, ...
                      'pass', e <= tolContact);

conv = [convShare(H5(sel, :)), convShare(LE(sel, :)), convShare(DM(sel, :)), convShare(cS(sel))]; % conv: change of the 3-sigma percentiles from N/2 to N
[e, je] = max(conv);
chk(end + 1) = struct('name', 'percentiles at half the draws against all, largest change over the 3-sigma spread', 'value', e, 'limit', tolConv, ...
                      'pass', e <= tolConv);

chk(end + 1) = struct('name', 'draw steps outside the drag table, height or LTAN held at its edge', 'value', sum(nOut), 'limit', 0, 'pass', sum(nOut) == 0);

R.checks = chk;
fprintf('checks, pass marks set before the first run:\n');
for j = 1:numel(chk)
    fprintf('   %2d. %s: %.4g, limit %.4g: %s\n', j, chk(j).name, chk(j).value, chk(j).limit, verdict(chk(j).pass));
end
fprintf('checks passed: %d of %d\n', nnz([chk.pass]), numel(chk));
R.nodeRuns = struct('ode', wRun, 'secular', dSec);
R.pathRuns = struct('ode', yOde(:, end).', 'cw', yCw(:, end).', 'dragOde', yDragOde(end), 'dragSecular', ySec);

% ---------------------------------------------------------------- figures for the slides
fig = cell(0, 3); % fig: figure name, width and height [in]
yb = tb/365.25; % yb: times of the stored percentiles [yr]
sz = struct('altitudeBand', [3.9 3.4], 'altitudeBandTop', [5.9 2.3], 'ltanBand', [5.9 2.3], 'shadowRange', [5.9 2.3], 'alongTrack', [5.9 2.3], ...
            'contactAltitude', [5.9 2.3], 'convergence', [4.4 2.2], 'decayCheck', [4.4 2.2], 'pathCheck', [3.0 3.0]); % sz: width and height of each figure [in]
hl = gobjects(1, 3); % hl: lines for the legends

% Altitude of the draws over the mission, median and 3-sigma band, each case labelled at its end:
nm = 'altitudeBand'; % nm: name of the figure being drawn
[f, ax] = slideFigure(nm, sz.(nm)(1), sz.(nm)(2), fs);
for i = 1:3
    fill(ax, [yb fliplr(yb)], [BH(1, :, i) fliplr(BH(3, :, i))], caseCol(i, :), 'FaceAlpha', 0.25, 'EdgeColor', 'none');
    plot(ax, yb, BH(2, :, i), 'Color', caseCol(i, :), 'LineWidth', 1.6);
end
plot(ax, [0 yrs], hLine*[1 1], 'k:', 'LineWidth', 1.0);
text(ax, 0.98*yrs, BH(3, end, 1) + 4, 'low', 'Color', caseCol(1, :), 'FontSize', fs, 'FontName', 'Arial', 'HorizontalAlignment', 'right', ...
     'VerticalAlignment', 'bottom');
text(ax, 0.98*yrs, BH(1, end, 2) - 4, 'moderate', 'Color', caseCol(2, :), 'FontSize', fs, 'FontName', 'Arial', 'HorizontalAlignment', 'right', ...
     'VerticalAlignment', 'top');
text(ax, 0.98*yrs, BH(2, end, 3), 'high', 'Color', caseCol(3, :), 'FontSize', fs, 'FontName', 'Arial', 'HorizontalAlignment', 'right', ...
     'VerticalAlignment', 'bottom');
set(ax, 'XLim', [0 yrs], 'XTick', 0:yrs, 'YLim', [base.decay.hFloor 10*ceil((max(BH(3, :, 1)) + 30)/10)]);
xlabel(ax, 'years in orbit');
ylabel(ax, 'orbit-mean altitude [km]');
exportgraphics(f, fullfile(figDir, [nm '.png']), 'Resolution', 300);
fig(end + 1, :) = {nm, sz.(nm)(1), sz.(nm)(2)};

% The same, ECSS low and moderate only:
nm = 'altitudeBandTop';
[f, ax] = slideFigure(nm, sz.(nm)(1), sz.(nm)(2), fs);
for i = 1:2
    fill(ax, [yb fliplr(yb)], [BH(1, :, i) fliplr(BH(3, :, i))], caseCol(i, :), 'FaceAlpha', 0.25, 'EdgeColor', 'none');
    plot(ax, yb, BH(2, :, i), 'Color', caseCol(i, :), 'LineWidth', 1.6);
end
plot(ax, [0 yrs], hLine*[1 1], 'k:', 'LineWidth', 1.0);
yLo = 10*floor(min([hLine, min(BH(1, :, 2))])/10 - 1); % yLo: bottom of the axis [km]
yHi = 10*ceil((max(BH(3, :, 1)) + 8)/10); % yHi: top of the axis, room for the label above the low band [km]
text(ax, 0.98*yrs, BH(3, end, 1) + 1, 'low', 'Color', caseCol(1, :), 'FontSize', fs, 'FontName', 'Arial', 'HorizontalAlignment', 'right', ...
     'VerticalAlignment', 'bottom');
text(ax, 0.98*yrs, BH(1, end, 2) - 1, 'moderate', 'Color', caseCol(2, :), 'FontSize', fs, 'FontName', 'Arial', 'HorizontalAlignment', 'right', ...
     'VerticalAlignment', 'top');
text(ax, 0.02*yrs, hLine - 1, sprintf('%g km', hLine), 'FontSize', fs, 'FontName', 'Arial', 'VerticalAlignment', 'top');
set(ax, 'XLim', [0 yrs], 'XTick', 0:yrs, 'YLim', [yLo yHi]);
xlabel(ax, 'years in orbit');
ylabel(ax, 'altitude [km]');
exportgraphics(f, fullfile(figDir, [nm '.png']), 'Resolution', 300);
fig(end + 1, :) = {nm, sz.(nm)(1), sz.(nm)(2)};

% LTAN of the draws over the mission, less the start:
nm = 'ltanBand';
[f, ax] = slideFigure(nm, sz.(nm)(1), sz.(nm)(2), fs);
for i = 1:3
    fill(ax, [yb fliplr(yb)], 60*[BL(1, :, i) fliplr(BL(3, :, i))], caseCol(i, :), 'FaceAlpha', 0.25, 'EdgeColor', 'none');
    plot(ax, yb, 60*BL(2, :, i), 'Color', caseCol(i, :), 'LineWidth', 1.6);
end
text(ax, 0.98*yrs, 60*BL(1, end, 1), 'low', 'Color', caseCol(1, :), 'FontSize', fs, 'FontName', 'Arial', 'HorizontalAlignment', 'right', ...
     'VerticalAlignment', 'top');
text(ax, 0.98*yrs, 60*BL(3, end, 2), 'moderate', 'Color', caseCol(2, :), 'FontSize', fs, 'FontName', 'Arial', 'HorizontalAlignment', 'right', ...
     'VerticalAlignment', 'bottom');
text(ax, 0.98*yrs, 60*BL(3, end, 3), 'high', 'Color', caseCol(3, :), 'FontSize', fs, 'FontName', 'Arial', 'HorizontalAlignment', 'right', ...
     'VerticalAlignment', 'top');
yr = 60*[min(reshape(BL(1, :, :), 1, [])), max(reshape(BL(3, :, :), 1, []))]; % yr: range of the bands [min]
set(ax, 'XLim', [0 yrs], 'XTick', 0:yrs, 'YLim', yr + [-0.15 0.08]*diff(yr));
xlabel(ax, 'years in orbit');
ylabel(ax, 'LTAN change [min]');
exportgraphics(f, fullfile(figDir, [nm '.png']), 'Resolution', 300);
fig(end + 1, :) = {nm, sz.(nm)(1), sz.(nm)(2)};

% Longest shadow over the mission, 3-sigma range of the draws in each case, median marked:
nm = 'shadowRange';
[f, ax] = slideFigure(nm, sz.(nm)(1), sz.(nm)(2), fs);
yy = [3 2 1]; % yy: row of each case
for i = 1:3
    S = R.levels(i).shadow; % S: nominal and percentiles [min]
    SW = R.levels(i).shadowWindow; % SW: the same with the 1-hour window
    plot(ax, [S(2) S(4)], yy(i) + 0.2*[1 1], '-', 'Color', caseCol(i, :), 'LineWidth', 7);
    plot(ax, S(3)*[1 1], yy(i) + 0.2 + [-0.12 0.12], 'k-', 'LineWidth', 1.5);
    plot(ax, [SW(2) SW(4)], yy(i) - 0.2*[1 1], '-', 'Color', 0.5*caseCol(i, :) + 0.5, 'LineWidth', 7);
    plot(ax, SW(3)*[1 1], yy(i) - 0.2 + [-0.12 0.12], 'k-', 'LineWidth', 1.5);
end
plot(ax, base.batt.tEcl*[1 1], [0.5 3.5], 'k:', 'LineWidth', 1.0);
S = R.levels(1).shadow; % S: top row again, for the labels left of its bars [min]
SW = R.levels(1).shadowWindow;
text(ax, S(2) - 0.4, yy(1) + 0.2, sprintf('RAAN %g^{\\circ}', mc.raan3), 'FontSize', fs, 'FontName', 'Arial', 'HorizontalAlignment', 'right', ...
     'VerticalAlignment', 'middle');
text(ax, SW(2) - 0.4, yy(1) - 0.2, sprintf('RAAN %g^{\\circ}', mc.raan3Window), 'FontSize', fs, 'FontName', 'Arial', 'HorizontalAlignment', 'right', ...
     'VerticalAlignment', 'middle');
xL = [5*floor(min(min(DM(sel, :)))/300 - 2), 5*ceil(max([base.batt.tEcl + 2, max(max(DW(sel, :)))/60])/5)]; % xL: axis range, room left of the bars for the labels [min]
set(ax, 'YLim', [0.5 3.5], 'YTick', 1:3, 'YTickLabel', fliplr(caseName), 'XLim', xL);
xlabel(ax, 'longest shadow [min]');
exportgraphics(f, fullfile(figDir, [nm '.png']), 'Resolution', 300);
fig(end + 1, :) = {nm, sz.(nm)(1), sz.(nm)(2)};

% Along-track error of a prediction against days ahead, 3 sigma:
td = 0:0.05:max(dAhead); % td: days ahead
gT = zeros(size(td)); % gT: GNSS fix [km]
dTk = zeros(3, numel(td)); % dTk: drag part, each case [km]
for a = 2:numel(td)
    gT(a) = pctl(abs(cw(td(a)*86400)), p3)/1e3;
    for i = 1:3
        dTk(i, a) = pctl(abs(yDrag(td(a)*86400, i)), p3)/1e3;
    end
end
nm = 'alongTrack';
[f, ax] = slideFigure(nm, sz.(nm)(1), sz.(nm)(2), fs);
plot(ax, td(2:end), gT(2:end), 'k-', 'LineWidth', 1.8);
for i = 1:3
    plot(ax, td(2:end), dTk(i, 2:end), '--', 'Color', caseCol(i, :), 'LineWidth', 1.6);
end
kg = round(0.35*(numel(td) - 1)) + 1; % kg: sample of the GNSS fix label, a third of the way along
text(ax, td(kg), gT(kg), 'GNSS fix', 'FontSize', fs, 'FontName', 'Arial', 'VerticalAlignment', 'bottom');
for i = 1:3
    text(ax, 0.98*max(dAhead), dTk(i, end), sprintf('drag, %s', caseName{i}), 'Color', caseCol(i, :), 'FontSize', fs, 'FontName', 'Arial', ...
         'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
end
set(ax, 'YScale', 'log', 'XLim', [0 max(dAhead)], 'XTick', 0:max(dAhead), 'YLim', [1e-2 10^ceil(log10(max(gT)))]);
xlabel(ax, 'days after the fix');
ylabel(ax, '3\sigma along-track [km]');
exportgraphics(f, fullfile(figDir, [nm '.png']), 'Resolution', 300);
fig(end + 1, :) = {nm, sz.(nm)(1), sz.(nm)(2)};

% Time in view a day against altitude:
nm = 'contactAltitude';
[f, ax] = slideFigure(nm, sz.(nm)(1), sz.(nm)(2), fs);
hh = linspace(min(hMean), max(hMean), 200); % hh: altitudes of the spline [km]
hl(1) = plot(ax, hh, cAt(hh), 'Color', col(1, :), 'LineWidth', 1.6);
hl(2) = plot(ax, hMean, cMin, 'o', 'Color', col(1, :), 'MarkerFaceColor', col(1, :), 'MarkerSize', 5);
hl(3) = plot(ax, [min(hMean) max(hMean)], comMin*[1 1], 'k:', 'LineWidth', 1.0);
xlabel(ax, 'orbit-mean altitude [km]');
ylabel(ax, 'in view [min/day]');
legend(ax, hl, {'spline', sprintf('ode45, %d days', nDayC), sprintf('%g min', comMin)}, 'Location', 'southeast', 'FontSize', fs);
exportgraphics(f, fullfile(figDir, [nm '.png']), 'Resolution', 300);
fig(end + 1, :) = {nm, sz.(nm)(1), sz.(nm)(2)};

% Convergence of the 3-sigma percentiles of the altitude after the mission with the number of draws:
nn = unique(round(logspace(3, log10(N), 14))); % nn: draw counts
nm = 'convergence';
[f, ax] = slideFigure(nm, sz.(nm)(1), sz.(nm)(2), fs);
for i = 1:3
    P = zeros(2, numel(nn)); % P: 0.135 and 99.865 percentiles with nn draws [km]
    for j = 1:numel(nn)
        P(:, j) = pctl(H5(1:nn(j), i), pct([1 3]));
    end
    sp = P(2, end) - P(1, end); % sp: 3-sigma spread with all draws [km]
    hl(i) = plot(ax, nn, 100*(P(1, :) - P(1, end))/sp, '-o', 'Color', caseCol(i, :), 'MarkerSize', 4, 'LineWidth', 1.3);
    plot(ax, nn, 100*(P(2, :) - P(2, end))/sp, '--s', 'Color', caseCol(i, :), 'MarkerSize', 4, 'LineWidth', 1.3);
end
plot(ax, [nn(1) nn(end)], 100*tolConv*[1 1], 'k:', [nn(1) nn(end)], -100*tolConv*[1 1], 'k:');
set(ax, 'XScale', 'log', 'XLim', [nn(1) nn(end)]);
xlabel(ax, 'draws');
ylabel(ax, 'change [%]');
legend(ax, hl, caseName, 'Location', 'northeast', 'FontSize', fs);
exportgraphics(f, fullfile(figDir, [nm '.png']), 'Resolution', 300);
fig(end + 1, :) = {nm, sz.(nm)(1), sz.(nm)(2)};

% Table decay of the sheet spacecraft less dragDecay.m:
nm = 'decayCheck';
[f, ax] = slideFigure(nm, sz.(nm)(1), sz.(nm)(2), fs);
for i = 1:3
    hd = M.dec{i, iS};
    m = min(numel(hd.h), size(hK, 2));
    hl(i) = plot(ax, (0:m-1)*dt/86400/365.25, hK(i, 1:m) - hd.h(1:m), 'Color', caseCol(i, :), 'LineWidth', 1.4);
end
set(ax, 'XLim', [0 yrs], 'XTick', 0:yrs);
xlabel(ax, 'years in orbit');
ylabel(ax, 'difference [km]');
legend(ax, hl, caseName, 'Location', 'southwest', 'FontSize', fs);
exportgraphics(f, fullfile(figDir, [nm '.png']), 'Resolution', 300);
fig(end + 1, :) = {nm, sz.(nm)(1), sz.(nm)(2)};

% Along-track after nDayP days, ode45 against Clohessy-Wiltshire:
nm = 'pathCheck';
[f, ax] = slideFigure(nm, sz.(nm)(1), sz.(nm)(2), fs);
lim = 1.1*max(abs([yOde(:, end); yCw(:, end)]))/1e3; % lim: axis half-range [km]
plot(ax, [-lim lim], [-lim lim], 'k:', 'LineWidth', 1.0);
plot(ax, yCw(:, end)/1e3, yOde(:, end)/1e3, 'o', 'Color', col(1, :), 'MarkerFaceColor', col(1, :), 'MarkerSize', 5);
set(ax, 'XLim', [-lim lim], 'YLim', [-lim lim]);
axis(ax, 'square');
xlabel(ax, 'Clohessy-Wiltshire [km]');
ylabel(ax, 'ode45 [km]');
exportgraphics(f, fullfile(figDir, [nm '.png']), 'Resolution', 300);
fig(end + 1, :) = {nm, sz.(nm)(1), sz.(nm)(2)};

R.figures = cell2struct(fig, {'name', 'w', 'h'}, 2).';
fprintf('\nfigures in %s: %s\n', figDir, strjoin(fig(:, 1).', ', '));

fid = fopen(fullfile(outDir, 'uncertainty.json'), 'w'); % fid: numbers file
fprintf(fid, '%s', jsonencode(R, 'PrettyPrint', true));
fclose(fid);
fprintf('numbers in %s and %s\n', txtFile, fullfile(outDir, 'uncertainty.json'));
diary off;


function tab = rateTable(cfg, hT, lT, tT, nUT, wSun)
% rateTable: orbit-averaged rate of change of semi-major axis at Cd A/m = 1 m^2/kg on a grid of height, LTAN offset and date, the orbit average of orbitRates in dragDecay.m, averaged over nUT start times spread over each day.
% cfg: parameter struct from config.m, activity level in cfg.drag
% hT: 1 x nH heights, semi-major axis less cfg.Re [km]
% lT: 1 x nL LTAN offsets, the node moved 15 deg/h east of one that follows the mean Sun from cfg.raan0 [h]
% tT: 1 x nT dates since sim start [days]
% nUT: start times a day
% wSun: mean Sun rate along the equator [deg/day]
% tab: struct of the grids h [km], L [h] and t [days], and lnK, nH x nL x nT natural log of minus the rate [log of km/s per m^2/kg]

if strcmp(cfg.drag.model, 'J77')
    cfg.drag.logRho = dragTable(cfg);
end

M = cfg.decay.nSample; % M: samples around each orbit
[H, L, U] = ndgrid(hT, lT, (0:nUT-1)/nUT); % H, L, U: height [km], LTAN offset [h] and start time of day [days] of each orbit of one date
a = cfg.Re + H(:).'; % a: 1 x B semi-major axes [km]
n = sqrt(cfg.mu./a.^3); % n: 1 x B mean motion [rad/s]
k2 = 1.5*cfg.J2*(cfg.Re./a).^2; % k2: 1 x B factor of the J2 secular rates
ci = cosd(cfg.inc); % ci: cosine of the inclination
si = sind(cfg.inc); % si: sine of the inclination

% Secular J2 rates of a circular orbit, node -k2 n cos i, argument of latitude n (1 + k2 (4 cos^2 i - 1)), as in dragDecay.m:
dW = -rad2deg(k2.*n*ci); % dW: 1 x B node rate [deg/s]
du = n.*(1 + k2*(4*ci^2 - 1)); % du: 1 x B rate of the argument of latitude [rad/s]
tau = (0:M-1).'*2*pi./du/M; % tau: M x B sample times after each start [s]

lnK = zeros(numel(hT), numel(lT), numel(tT));
for k = 1:numel(tT)
    t0 = 86400*(tT(k) + U(:).'); % t0: 1 x B start times [s]
    W = cfg.raan0 + wSun*t0/86400 + 15*L(:).'; % W: 1 x B node at each start [deg]
    u = deg2rad(cfg.argp0 + cfg.nu0) + du.*(t0 + tau); % u: M x B argument of latitude, carried at the J2 rate [rad]
    Om = deg2rad(W + dW.*tau); % Om: M x B node [rad]

    % Circular orbit, position and velocity from the node, inclination and argument of latitude, as in dragDecay.m:
    r = [a.*(cos(Om).*cos(u) - sin(Om).*sin(u)*ci); a.*(sin(Om).*cos(u) + cos(Om).*sin(u)*ci); a.*sin(u)*si];
    % r: 3M x B positions, x rows 1 to M, y rows M+1 to 2M, z rows 2M+1 to 3M [km]
    vc = sqrt(cfg.mu./a); % vc: 1 x B circular speed [km/s]
    v = [vc.*(-cos(Om).*sin(u) - sin(Om).*cos(u)*ci); vc.*(-sin(Om).*sin(u) + cos(Om).*cos(u)*ci); vc.*cos(u)*si];
    % v: 3M x B velocities, rows as r [km/s]
    r3 = [reshape(r(1:M, :), 1, []); reshape(r(M+1:2*M, :), 1, []); reshape(r(2*M+1:end, :), 1, [])]; % r3: 3 x MB positions [km]
    v3 = [reshape(v(1:M, :), 1, []); reshape(v(M+1:2*M, :), 1, []); reshape(v(2*M+1:end, :), 1, [])]; % v3: 3 x MB velocities [km/s]

    rho = dragDensity(reshape(t0 + tau, 1, []), r3, cfg); % rho: 1 x MB mass density [kg/m^3]

    % Velocity relative to an atmosphere rotating rigidly with the Earth, v - wE x r:
    vRel = v3 - cfg.wE*[-r3(2, :); r3(1, :); zeros(1, size(r3, 2))]; % vRel: 3 x MB [km/s]

    % Drag acceleration at Cd A/m = 1 m^2/kg, 1e3 converting (kg/m^3)(m^2/kg)(km/s)^2 to km/s^2:
    aD = -0.5e3*rho.*vecnorm(vRel).*vRel; % aD: 3 x MB [km/s^2]

    % Semi-major axis rate from the specific energy rate, da/dt = (2 a^2/mu) <v . aD>:
    dadt = 2*a.^2/cfg.mu.*mean(reshape(sum(v3.*aD, 1), M, []), 1); % dadt: 1 x B [km/s]

    lnK(:, :, k) = log(-mean(reshape(dadt, numel(hT), numel(lT), nUT), 3));
end

tab.h = hT(:);
tab.L = lT(:).';
tab.t = tT(:).';
tab.lnK = lnK;

end


function [K, out] = tableRate(tab, h, L, t)
% tableRate: orbit-averaged rate of change of semi-major axis at Cd A/m = 1 m^2/kg from the drag table, linear in date and, in the log of the rate, in height and LTAN offset.
% tab: drag table struct from rateTable
% h: N x 1 heights [km], held at the edges of the table
% L: N x 1 LTAN offsets [h], held at the edges of the table
% t: time since sim start [s]
% K: N x 1 rate [km/s per m^2/kg]
% out: N x 1 true where h or L was held at an edge

x = t/86400/(tab.t(2) - tab.t(1)); % x: date in table steps
j = min(floor(x) + 1, numel(tab.t) - 1); % j: date before t
w = x - (j - 1); % w: share of the way to the next date
S = (1 - w)*tab.lnK(:, :, j) + w*tab.lnK(:, :, j + 1); % S: nH x nL log rate at t
Lc = min(max(L, tab.L(1)), tab.L(end)); % Lc: LTAN offset held inside the table [h]
hc = min(max(h, tab.h(1)), tab.h(end)); % hc: height held inside the table [km]
out = Lc ~= L | hc ~= h;
K = -exp(interp2(tab.L, tab.h, S, Lc, hc, 'linear'));

end


function [hS, L, tF, nOut] = decayStep(hS, L, tF, t, F, inc, tab, rateL, cfg)
% decayStep: one midpoint step of dragDecay.m for many circular orbits, the orbit-averaged rate from the drag table, the LTAN offset carried along.
% hS: N x 1 altitude, semi-major axis less cfg.Re in the circular model [km]
% L: N x 1 LTAN offset [h]
% tF: N x 1 time each orbit fell below cfg.decay.hFloor, NaN while above [s]
% t: time at the start of the step [s]
% F: N x 1 Cd A/m times the density factor [m^2/kg]
% inc: N x 1 inclination [deg]
% tab: drag table struct from rateTable
% rateL: handle giving the LTAN rate [h/day] from the altitude [km] and inclination [deg]
% cfg: parameter struct from config.m
% nOut: orbits whose height or LTAN offset was held at an edge of the drag table in this step

dt = 86400*cfg.decay.dt; % dt: step length [s]
q = isnan(tF); % q: orbits still above the floor
h = hS(q); % h: their altitude [km]
l = L(q); % l: their LTAN offset [h]
f = F(q); % f: their drag factor [m^2/kg]
in = inc(q); % in: their inclination [deg]

[K1, o1] = tableRate(tab, h, l, t); % K1: rate at the start of the step [km/s per m^2/kg]
dh1 = f.*K1; % dh1: rate of change of altitude [km/s]
dl1 = rateL(h, in)/86400; % dl1: rate of change of the LTAN offset [h/s]

% Midpoint rule, as in dragDecay.m:
[K2, o2] = tableRate(tab, h + dt/2*dh1, l + dt/2*dl1, t + dt/2);
hN = h + dt*f.*K2; % hN: altitude at the end of the step [km]
lN = l + dt*rateL(h + dt/2*dh1, in)/86400; % lN: LTAN offset at the end of the step [h]

% Floor crossing, linear between the two ends of the step as in dragDecay.m:
fell = hN < cfg.decay.hFloor; % fell: orbits below the floor at the end of the step
tq = tF(q); % tq: floor times of the orbits still above it [s]
tq(fell) = t + (cfg.decay.hFloor - h(fell))./(hN(fell) - h(fell))*dt;
tF(q) = tq;
hS(q) = hN;
L(q) = lN;
nOut = nnz(o1 | o2);

end


function phi = shadowTable(cfg, bP, hP)
% shadowTable: half-width of the Earth shadow arc of a circular orbit about the point opposite the Sun, against the Sun angle above the orbit plane and the height, penumbra counted as in eclipseTimes.m, Sun at 1 AU, by bisection on the penumbra margin of shadowGeometry.m.
% cfg: parameter struct from config.m
% bP: Sun angles [deg]
% hP: heights, orbit radius less cfg.Re [km]
% phi: struct of the grids b [deg] and h [km], and half, numel(bP) x numel(hP) half-width, zero without shadow [rad]

[B, H] = ndgrid(bP, hP); % B: Sun angle [deg], H: height [km] of each point
r0 = cfg.Re + H(:).'; % r0: 1 x P orbit radius [km]

% Sun in the frame of the orbit, x toward its projection on the orbit plane, z along the orbit normal:
rs = cfg.AU*[cosd(B(:).'); zeros(1, numel(B)); sind(B(:).')]; % rs: 3 x P Sun position [km]

lo = zeros(1, numel(B)); % lo: half-widths still in shadow [rad]
hi = pi*ones(1, numel(B)); % hi: half-widths in sunlight [rad]
in = penumbra(pi*ones(1, numel(B)), r0, rs, cfg) < 0; % in: points with a shadow
for j = 1:45
    mid = (lo + hi)/2;
    q = penumbra(pi - mid, r0, rs, cfg) < 0; % q: still in shadow mid away from the point opposite the Sun
    lo(q) = mid(q);
    hi(~q) = mid(~q);
end
half = (lo + hi)/2; % half: half-width [rad]
half(~in) = 0;

phi.b = bP(:);
phi.h = hP(:).';
phi.half = reshape(half, numel(bP), numel(hP));

end


function g = penumbra(u, r0, rs, cfg)
% penumbra: penumbra margin of eclipseTimes.m at an argument of latitude measured from the Sun's projection on the orbit plane.
% u: 1 x P argument of latitude [rad]
% r0: 1 x P orbit radius [km]
% rs: 3 x P Sun position in the orbit frame [km]
% cfg: parameter struct from config.m
% g: 1 x P margin, negative inside the penumbra [rad]

r = [r0.*cos(u); r0.*sin(u); zeros(size(u))]; % r: 3 x P satellite position [km]
[sep, aSun, aEarth] = shadowGeometry(r, rs, cfg);
g = sep - (aEarth + aSun);

end


function [D, beta] = shadowLength(s, h, inc, W, uDot, phi)
% shadowLength: Earth shadow of circular orbits, penumbra entry to exit, from the shadow table at their Sun angle and height.
% s: 3 x 1 unit vector to the Sun
% h: N x 1 orbit-mean altitude [km]
% inc: N x 1 inclination [deg]
% W: N x 1 node [deg]
% uDot: N x 1 rate of the argument of latitude [rad/s]
% phi: shadow table struct from shadowTable
% D: N x 1 shadow length [s]
% beta: N x 1 Sun angle above the orbit plane [deg]

% Orbit normal from the node and inclination:
hn = [sind(W).*sind(inc), -cosd(W).*sind(inc), cosd(inc).*ones(size(W))]; % hn: N x 3 unit orbit normal

% Sun elevation above the orbit plane, beta = asin(h_hat . s_hat):
beta = asind(hn*s(:));

D = 2*interp2(phi.h, phi.b, phi.half, h, abs(beta), 'linear')./uDot;

end


function [t, X] = odeRun(cfg, x0, nDay)
% odeRun: the integration of propagate.m from a given state, output every cfg.ode.dtOut.
% cfg: parameter struct from config.m
% x0: 6 x 1 initial state [km, km/s]
% nDay: run length [days]
% t: output times [s]
% X: states, one row per output time [km, km/s]

if cfg.use.drag
    cfg.drag.logRho = dragTable(cfg);
end
opts = odeset('RelTol', cfg.ode.relTol, 'AbsTol', cfg.ode.absTol); % opts: ode45 options structure
[t, X] = ode45(@(t, x) eom(t, x, cfg), 0:cfg.ode.dtOut:nDay*86400, x0, opts);

end


function y = alongTrack(X0, X1)
% alongTrack: along-track distance of one orbit ahead of another, the angle between the positions about the first orbit's normal times its radius.
% X0: K x 6 states of the reference orbit [km, km/s]
% X1: K x 6 states of the other orbit at the same times [km, km/s]
% y: 1 x K along-track distance, positive ahead [m]

hn = cross(X0(:, 1:3), X0(:, 4:6), 2); % hn: K x 3 angular momentum of the reference [km^2/s]
hn = hn./vecnorm(hn, 2, 2);
ang = atan2(sum(cross(X0(:, 1:3), X1(:, 1:3), 2).*hn, 2), sum(X0(:, 1:3).*X1(:, 1:3), 2)); % ang: K x 1 angle from the reference [rad]
y = 1e3*(vecnorm(X0(:, 1:3), 2, 2).*ang).';

end


function [cm, hm, cd] = contactRun(cfg, nDay, st)
% contactRun: time in view of any station a day over a run of nDay days, from downlinkTimes.m.
% cfg: parameter struct from config.m
% nDay: run length [days]
% st: station structs from stations.m
% cm: mean time in view a day [min]
% hm: orbit-mean altitude over the run [km]
% cd: 1 x nDay time in view each day [min]

% Kepler's third law:
T = 2*pi*sqrt((cfg.Re + cfg.h)^3/cfg.mu); % T: orbital period [s]

out = propagate(cfg, (nDay*86400 + 1)/T); % out: state history over nDay days
dl = downlinkTimes(out, cfg, st); % dl: windows of any station, start and end [s]
cd = zeros(1, nDay);
for d = 1:nDay
    cd(d) = sum(max(0, min(dl(:, 2), d*86400) - max(dl(:, 1), (d - 1)*86400)))/60;
end
cm = mean(cd);
hm = mean(vecnorm(out.r)) - cfg.Re;

end


function p = pctl(x, q)
% pctl: percentiles of a sample, linear between the order statistics placed at (k - 0.5)/n as prctile does, Inf kept.
% x: sample
% q: percentiles [%]
% p: 1 x numel(q) values

x = sort(x(:)); % x: sorted sample
n = numel(x); % n: sample size
pos = min(max(q(:).'/100*n + 0.5, 1), n); % pos: position of each percentile in the sorted sample
lo = floor(pos); % lo: order statistic below
w = pos - lo; % w: share of the way to the one above
p = x(lo).';
k = w > 0; % k: percentiles between two order statistics
p(k) = x(lo(k)).'.*(1 - w(k)) + x(lo(k) + 1).'.*w(k);

end


function e = convShare(X)
% convShare: largest change of the 3-sigma percentiles from the first half of the draws to all, over the 3-sigma spread, each column.
% X: N x C draws
% e: 1 x C change, share of the spread

n = size(X, 1); % n: draws
e = zeros(1, size(X, 2));
for j = 1:size(X, 2)
    P1 = pctl(X(1:floor(n/2), j), [0.135 99.865]); % P1: percentiles of the first half
    P2 = pctl(X(:, j), [0.135 99.865]); % P2: percentiles of all
    sp = P2(2) - P2(1); % sp: 3-sigma spread
    if sp > 0 && all(isfinite([P1 P2]))
        e(j) = max(abs(P1 - P2))/sp;
    end
end

end


function s = verdict(tf)
% verdict: PASS or FAIL for a check.
% tf: check result
% s: text for the log

if tf
    s = 'PASS';
else
    s = 'FAIL';
end

end
