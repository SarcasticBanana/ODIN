% MISSION: the five-year mission from the Inputs sheet, everything on the DR2 power and drag slides: shadow of every orbit of each plane, the worst orbit, array power, battery ageing, energy margin and lowest charge of each operating mode at end of life, observation windows, density by local time and drag decay against m/(Cd A). Figures in results/mission at 14 pt, every number in results/mission.txt and results/mission.json, the runs in results/mission.mat.

clear; close all; clc;

here = fileparts(mfilename('fullpath')); % here: repository root folder
addpath(here);
odinPath;

rerun = true; % rerun: false redraws the figures and rewrites the numbers from results/mission.mat without the runs
slewMin = [30 45 60]; % slewMin: observation window lengths [min]
tiltMax = 90; % tiltMax: worst-case window tilt from the Sun, ADCS [deg]
Bg = [10 15 20 30 40 60 80 100 150 200 300 400 600 1000]; % Bg: ballistic coefficients m/(Cd A) of the plane 1 decay runs [kg/m^2]
Bp = [30 100]; % Bp: ballistic coefficients of the plane 2 and 3 decay runs [kg/m^2]
hLine = 550; % hLine: altitude the five-year decay is judged against, the lowest the battery is sized for (cfg.batt.h) [km]
fHz = 1; % fHz: EPS telemetry frame rate [Hz]
fs = 14; % fs: font size of every figure, the 14 pt floor of the DR2 assignment [pt]
frame = {'time stamp (32-bit) and status flags', 2, 4
         'load channels V and I (5 channels)', 10, 2
         'array inputs V and I (2 wings)', 4, 2
         'battery V, I and 2 temperatures (2 packs)', 8, 2
         'cell voltages (2 x 4S) and charge', 9, 2
         'bus voltage, board temperature', 2, 2}; % frame: EPS telemetry frame from the EPS block diagram channels, item, values and bytes per value

outDir = fullfile(here, 'results'); % outDir: folder for the numbers and the saved runs
figDir = fullfile(outDir, 'mission'); % figDir: folder for the slide figures
if ~isfolder(figDir)
    mkdir(figDir);
end
matFile = fullfile(outDir, 'mission.mat'); % matFile: saved runs

col = [0.00 0.45 0.70
       0.84 0.37 0.00
       0.00 0.62 0.45
       0.80 0.47 0.65
       0.90 0.62 0.00
       0.34 0.71 0.91]; % col: 6 x 3 plane colours, as in RUN.m
caseCol = col([3 5 2], :); % caseCol: 3 x 3 colour of the ECSS low, moderate and high cases, as in RUN.m
caseName = {'low', 'moderate', 'high'}; % caseName: ECSS activity level of each case
modeName = {'Nominal', 'Raven', 'Safe'}; % modeName: operating modes
grey = [0.5 0.5 0.5]; % grey: colour of reference lines

base = config(); % base: parameter struct of ODIN 1.1 from the Inputs sheet
[cfgs, names] = constellation(base);
% cfgs: parameter struct of each satellite from constellation.m
% names: name of each satellite

nP = base.const.nPerPlane; % nP: satellites per plane
nPl = base.const.nPlanes; % nPl: number of orbit planes
yrs = base.power.years; % yrs: mission length [yr]
nDay = round(365.25*yrs); % nDay: days in the mission
dn0 = datenum(base.epoch); % dn0: date number at sim start [days]
F107 = [base.decay.F107Low, base.decay.F107Moderate, base.decay.F107High]; % F107: 10.7 cm solar flux of the ECSS low, moderate and high cases [sfu]
Ap = [base.decay.ApLow, base.decay.ApModerate, base.decay.ApHigh]; % Ap: daily geomagnetic index of the same cases
pL = [base.modes.pNominal, base.modes.pRaven, base.modes.pSafe]; % pL: bus load of each mode [W]
tilts = [base.slew.tilt, tiltMax]; % tilts: window tilts from the Sun [deg]
Bs = base.drag.mass/(base.drag.Cd*base.drag.area); % Bs: ballistic coefficient of the spacecraft in the Inputs sheet [kg/m^2]
Ball = [Bg, Bs, 1e6]; % Ball: ballistic coefficients of the plane 1 runs, then the sheet spacecraft, the last so large the altitude holds [kg/m^2]
dateOf = @(t) datestr(dn0 + t/86400, 'dd mmm yyyy'); % dateOf: date text of a time since sim start

if rerun
    tic;

    sv = cell(1, nPl); % sv: shadowSurvey of satellite k.1 of each plane, drag off so the altitude holds
    for k = 1:nPl
        c = cfgs((k - 1)*nP + 1); % c: parameter struct of satellite k.1
        c.use.drag = false;
        sv{k} = shadowSurvey(c, nDay);
        fprintf('plane %d: %d shadows in %d days, longest %.2f min, %.0f s so far\n', k, numel(sv{k}.dur), nDay, sv{k}.dMax/60, toc);
    end
    [~, kW] = max(cellfun(@(s) s.dMax, sv)); % kW: plane with the longest shadow
    if sv{kW}.dMax == 0
        error('MISSION:noShadow', 'no Earth shadow in %g years, the battery run needs one', yrs);
    end

    % Battery through the longest shadow every orbit at the Nominal load, the whole mission:
    cB = base; % cB: parameter struct of the battery run
    cB.batt.pEcl = base.modes.pNominal;
    cB.batt.tEcl = sv{kW}.dMax/60;
    cB.batt.h = sv{kW}.h0;
    b = batteryAging(cB); % b: battery capacity, resistance and lowest charge from batteryAging.m

    % End-of-life battery, the capacity and discharge efficiency left after the mission:
    cE = base; % cE: parameter struct with the end-of-life battery
    cE.batt.energy = base.batt.energy*b.qEOL;
    cE.batt.etaDis = b.etaDisEOL;

    S = sv{kW}.first; % S: state history over the first week
    W = sv{kW}.worst; % W: state history over the worst week
    pmS = {powerModes(S, cE, base.slew.tiltSun), powerModes(S, cE, base.slew.tilt)}; % pmS: first week on the Sun, and tilted by cfg.slew.tilt all week
    pmW = {powerModes(W, cE, base.slew.tiltSun), powerModes(W, cE, base.slew.tilt)}; % pmW: the same over the worst week
    frS = slewFrequency(S, cE, slewMin, tilts); % frS: observation windows in the first week
    frW = slewFrequency(W, cE, slewMin, tilts); % frW: observation windows in the worst week
    fprintf('battery, power modes and observation windows done, %.0f s so far\n', toc);

    % Density around the equator at the orbit-mean altitude by local solar time at sim start, as in the Drag day and night figure of RUN.m:
    th0 = earthRotation(0, base); % th0: Greenwich angle at sim start [rad]
    ut0 = 3600*base.epoch(4) + 60*base.epoch(5) + base.epoch(6); % ut0: UTC time of day at sim start [s]
    lstG = 0:0.25:24; % lstG: local solar times around the equator [h]
    lonG = deg2rad(15*lstG - ut0/240); % lonG: longitudes where the local solar time is lstG at sim start [rad]
    rG = (base.Re + sv{1}.h0)*[cos(lonG + th0); sin(lonG + th0); zeros(size(lonG))]; % rG: 3 x L equator points in the model frame [km]
    rhoG = zeros(3, numel(lstG)); % rhoG: 3 x L NRLMSISE-00 density in each case [kg/m^3]
    for i = 1:3
        cN = base; % cN: parameter struct of case i with NRLMSISE-00
        cN.drag.model = 'NRLMSISE-00';
        cN.drag.F107 = F107(i);
        cN.drag.Ap = Ap(i);
        rhoG(i, :) = dragDensity(zeros(size(lstG)), rG, cN);
    end

    % Decay of plane 1 from its orbit-mean altitude in each case at each ballistic coefficient, B set through the mass with Cd A = 1 m^2:
    dec = cell(3, numel(Ball)); % dec: dragDecay of satellite 1.1 in case i at Ball(j)
    for i = 1:3
        for j = 1:numel(Ball)
            c = cfgs(1);
            c.h = sv{1}.h0;
            c.drag.F107 = F107(i);
            c.drag.Ap = Ap(i);
            c.drag.mass = Ball(j);
            c.drag.Cd = 1;
            c.drag.area = 1;
            dec{i, j} = dragDecay(c, 365.25*yrs);
        end
        fprintf('decay at ECSS %s done, %.0f s so far\n', caseName{i}, toc);
    end
    decP = cell(3, nPl, numel(Bp)); % decP: dragDecay of satellite k.1 in case i at Bp(j), planes 2 onward
    for i = 1:3
        for k = 2:nPl
            for j = 1:numel(Bp)
                c = cfgs((k - 1)*nP + 1);
                c.h = sv{k}.h0;
                c.drag.F107 = F107(i);
                c.drag.Ap = Ap(i);
                c.drag.mass = Bp(j);
                c.drag.Cd = 1;
                c.drag.area = 1;
                decP{i, k, j} = dragDecay(c, 365.25*yrs);
            end
        end
    end
    fprintf('all runs done in %.0f s\n', toc);

    save(matFile, 'sv', 'kW', 'b', 'cB', 'cE', 'pmS', 'pmW', 'frS', 'frW', 'lstG', 'rhoG', 'dec', 'decP', '-v7');
else
    load(matFile);
end

S = sv{kW}.first; % S: state history over the first week
W = sv{kW}.worst; % W: state history over the worst week
TN = sv{kW}.TN; % TN: nodal period [s]
tEcl = sv{kW}.dMax; % tEcl: longest shadow over the mission [s]

% ---------------------------------------------------------------- orbit and shadows
wSun = 36000.771/36525; % wSun: mean Sun rate along the equator, as in the LTAN of RUN.m [deg/day]
R.orbit.h = base.h; % osculating altitude at sim start [km]
R.orbit.inc = base.inc; % inclination [deg]
R.orbit.h0 = sv{1}.h0; % orbit-mean altitude [km]
R.orbit.TN = TN; % nodal period [s]
R.orbit.dW = sv{1}.dW; % node rate of the ode45 orbit [deg/day]
R.orbit.wSun = wSun;
R.orbit.ltanDrift = (sv{1}.dW - wSun)*365.25*4; % LTAN drift, 4 min per degree [min/yr]

% Sun-synchronous inclination, node rate proportional to cos i as in ssoInclination of tradeLighting.m:
R.orbit.incSSO = acosd(cosd(base.inc)*(360/365.2422)/sv{1}.dW);

R.orbit.ltan0 = sv{1}.ltan(1); % LTAN of plane 1 at sim start [h]
R.orbit.raan0 = base.raan0;
R.orbit.dRaan = base.const.dRaan;
R.orbit.epoch = datestr(dn0, 'dd mmm yyyy');
R.years = yrs;
R.nDay = nDay;
R.worstPlane = kW;
for k = 1:nPl
    P = struct(); % P: shadow numbers of plane k
    P.longest = sv{k}.dMax/60; % longest shadow [min]
    P.date = dateOf(sv{k}.tMax);
    P.tMax = sv{k}.tMax; % entry of the longest shadow [s]
    P.beta = sv{k}.betaMax; % beta angle there [deg]
    P.ltan = sv{k}.ltanMax; % LTAN there [h]
    P.shadowFree = nnz(sv{k}.day == 0); % days without shadow
    P.firstShadowDay = find(sv{k}.day > 0, 1) - 1; % first day with a shadow, 0 on sim start
    P.h0 = sv{k}.h0;
    P.ltan0 = sv{k}.ltan(1); % LTAN at sim start [h]
    P.betaMin = min(abs(sv{k}.beta)); % smallest Sun angle above or below the orbit plane at the day starts [deg]
    P.betaMax = max(abs(sv{k}.beta)); % largest [deg]
    P.yearLongest = zeros(1, ceil(yrs)); % longest shadow of each mission year [min]
    P.yearShare = zeros(1, ceil(yrs)); % share of the days of each mission year with a shadow
    for y = 1:ceil(yrs)
        qd = (1:nDay) > (y - 1)*365.25 & (1:nDay) <= y*365.25; % qd: days of year y
        P.yearLongest(y) = max(sv{k}.day(qd))/60;
        P.yearShare(y) = nnz(sv{k}.day(qd) > 0)/nnz(qd);
    end
    R.planes(k) = P;
end
R.worst.date = dateOf(sv{kW}.tMax);
R.worst.dur = tEcl/60; % longest shadow [min]
R.worst.beta = sv{kW}.betaMax;
R.worst.ltan = sv{kW}.ltanMax;
R.worst.weekStart = dateOf(W.t(1));
R.worst.weekEnd = dateOf(W.t(end));
R.first.weekStart = dateOf(S.t(1));
R.first.weekEnd = dateOf(S.t(end));
R.battSizing.h = base.batt.h; % battery sizing orbit of the sheet [km]
R.battSizing.tEcl = base.batt.tEcl; % its longest shadow [min]

% ---------------------------------------------------------------- loads, currents, limits and telemetry
ld = base.load; % ld: component loads from the Loads section of the sheet [W]
R.loads.p = pL; % bus load of each mode [W]
R.loads.pRelease = ld.release; % burn wires at array release [W]
R.loads.component = {'C&DH board', 'ADCS (XACT-50)', 'Event camera', 'Radio (AL-TOTEM)'};
R.loads.pComponent = [ld.cdh, ld.adcs, ld.camera, ld.radio]; % component loads at the top of each range [W]
R.loads.iComponent = R.loads.pComponent/base.batt.vMin; % component currents at the lowest battery voltage [A]
R.loads.iNominal = base.modes.pNominal/base.batt.vMin;
R.loads.vMin = base.batt.vMin;
R.loads.sumCheck = round(1e6*([ld.cdh + ld.adcs + ld.camera + ld.radio, ld.cdh + ld.adcs + ld.radio, ld.cdh + ld.adcs + ld.radio] - pL))/1e6;
% R.loads.sumCheck: components on in each mode less the mode load, zero when the sheet is consistent [W]
R.limits.iChannel = base.eps.iChannel; % PDU channel limit [A]
R.limits.iPack = base.batt.iPack; % battery pack discharge limit [A]
R.limits.nPack = base.batt.nPack;
R.limits.packMass = base.batt.packMass; % mass of one pack [kg]
R.limits.pBatt = base.batt.nPack*base.batt.iPack*base.batt.vMin; % battery power limit at the lowest voltage [W]
R.limits.pPeak = [pL(1), pL(2) + ld.release, pL(3)]; % peak load of each mode [W]
R.limits.heater = base.batt.vMax^2/base.batt.heaterR; % one pack heater at the top of charge [W]
R.limits.epsIdle = base.eps.pIdle; % EPS electronics idle [W]
R.frame.item = frame(:, 1).';
R.frame.values = cell2mat(frame(:, 2)).';
R.frame.bytes = cell2mat(frame(:, 2)).'.*cell2mat(frame(:, 3)).';
R.frame.total = sum(R.frame.bytes); % bytes per frame
R.frame.bps = 8*R.frame.total*fHz; % bit rate [bit/s]
R.frame.MBday = R.frame.total*fHz*86400/1e6; % data a day [MB]
R.frame.fHz = fHz;

% ---------------------------------------------------------------- battery
eFull0 = base.batt.energy*base.batt.socTop; % eFull0: full charge at the start of life [Wh]
R.batt.energy = base.batt.energy;
R.batt.socTop = base.batt.socTop;
R.batt.vMin = base.batt.vMin;
R.batt.vMax = base.batt.vMax;
R.batt.vNom = base.batt.vNom;
R.batt.eFull0 = eFull0;
R.batt.eFullEOL = eFull0*b.qEOL; % full charge at the end of life [Wh]
R.batt.eDraw = b.eDraw; % battery energy through the longest shadow at the start of life [Wh]
R.batt.kR = b.kR;
R.batt.qEOL = b.qEOL;
R.batt.rEOL = b.rEOL;
R.batt.etaDis = base.batt.etaDis;
R.batt.etaDisEOL = b.etaDisEOL;
R.batt.etaChg = base.batt.etaChg;
R.batt.efc = b.efc(end); % equivalent full cycles over the mission
R.batt.socMinBOL = b.socMin(1);
R.batt.socMinEOL = b.socMin(end);
R.batt.tFail = b.tFail;
R.batt.floor = base.batt.floor;
R.batt.temp = base.batt.temp;
R.batt.fadeScale = base.batt.fadeScale;
R.batt.dose = base.batt.dose;
R.batt.model = base.batt.model;
for y = 0:ceil(yrs)
    kq = min(round(365.25*y) + 1, numel(b.t)); % kq: day index at the start of year y
    R.batt.yearQ(y + 1) = b.q(kq);
    R.batt.yearR(y + 1) = b.r(kq);
    R.batt.yearSoc(y + 1) = b.socMin(kq);
end

% Lowest charge of each mode through the longest shadow, 1 - P t/(eta_d E_full):
R.batt.socModeBOL = 1 - pL*tEcl/3600/base.batt.etaDis/eFull0;
R.batt.socModeEOL = 1 - pL*tEcl/3600/b.etaDisEOL/(eFull0*b.qEOL);

% ---------------------------------------------------------------- array over one orbit, the last nodal period of each week
r = arrayPower(struct('t', W.t(1:2), 'illum', W.illum(1:2), 'rsun', W.rsun(:, 1:2)), base, base.power.panelArea, ...
               struct('incidenceAngles_deg', 0, 'incidenceSummaryAngles_deg', 0, 'warnTemperatureExtrapolation', false, 'warnLumpedModel', false));
% r: arrayPower over one step, run for its parameter struct
p = r.params; % p: arrayPower parameters
caseArr = {'first week, on the Sun', sprintf('first week, %g deg', base.slew.tilt), 'worst orbit, on the Sun', sprintf('worst orbit, %g deg', base.slew.tilt)};
% caseArr: label of each array case
oArr = {S, S, W, W}; % oArr: state history of each case
pmArr = {pmS{1}, pmS{2}, pmW{1}, pmW{2}}; % pmArr: powerModes run of each case
A = struct('t', {}, 'p', {}, 'T', {}); % A: one orbit of each case
for c = 1:4
    o = oArr{c}; % o: state history of case c
    tA = o.t(end) - TN; % tA: start of the last nodal period [s]
    q = o.t > tA; % q: samples after tA
    A(c).t = [tA, o.t(q)] - tA; % time into the orbit [s]
    A(c).p = [interp1(o.t, pmArr{c}.pArr, tA), pmArr{c}.pArr(q)]; % array power at the bus [W]
    A(c).T = [interp1(o.t, pmArr{c}.T, tA), pmArr{c}.T(q)] - 273.15; % panel temperature [deg C]
    A(c).eGen = trapz(A(c).t, A(c).p)/3600; % energy over the orbit [Wh]
    A(c).mean = A(c).eGen*3600/TN; % mean power [W]
    A(c).peak = max(A(c).p); % peak power [W]
    A(c).Tmin = min(A(c).T);
    A(c).Tmax = max(A(c).T);

    % Temperature coefficient held inside its fit range, power scaled by (1 + g (Tc - Tref))/(1 + g (T - Tref)):
    Tc = min(max(A(c).T, p.powerTempCoeffRange_C(1)), p.powerTempCoeffRange_C(2)); % Tc: temperature held inside the range [deg C]
    fac = (1 + p.powerTempCoeff_perC*(Tc - p.referenceTemperature_C))./(1 + p.powerTempCoeff_perC*(A(c).T - p.referenceTemperature_C));
    % fac: power with the held temperature over the model power
    A(c).eClamp = trapz(A(c).t, A(c).p.*fac)/3600; % energy over the orbit with the held temperature [Wh]
    A(c).tSh = sum(sv{kW}.dur(sv{kW}.in >= o.t(end) - TN & sv{kW}.in < o.t(end))); % shadow inside the orbit [s]
    A(c).name = caseArr{c};
end
R.array.panelArea = base.power.panelArea;
R.array.eol = p.agingRetained;
R.array.cellEff = p.cellEfficiency;
R.array.tempCoeff = p.powerTempCoeff_perC;
R.array.tempRef = p.referenceTemperature_C;
R.array.tempRange = p.powerTempCoeffRange_C;
R.array.tempStart = p.initialTemperature_C;
R.array.cases = rmfield(A, {'t', 'p', 'T'});

% ---------------------------------------------------------------- energy margin of each mode, Sun-pointed, end-of-life battery
mS = energyMargin(A(1).t, A(1).p, pL, A(1).tSh, cE); % mS: first week
mW = energyMargin(A(3).t, A(3).p, pL, A(3).tSh, cE); % mW: worst orbit
R.margin.first = mS;
R.margin.worst = mW;
R.week.socMin = min(pmW{1}.soc); % lowest charge over the worst week, Nominal on the Sun
R.week.hMode = pmW{1}.tMode/3600; % hours in Nominal, Raven and Safe

% ---------------------------------------------------------------- observation windows
R.slew.first = struct('tilt', {}, 'D', {}, 'perDay', {}, 'recMed', {}, 'recMax', {}, 'socMin', {}, 'hMode', {}, 'nCut', {});
R.slew.worst = R.slew.first;
for i = 1:numel(tilts)
    for j = 1:numel(slewMin)
        for o = 1:2
            if o == 1
                f = frS(i, j); % f: one case of slewFrequency
            else
                f = frW(i, j);
            end
            e = struct('tilt', f.tilt, 'D', f.D, 'perDay', f.perDay, 'recMed', median(f.rec(~isnan(f.rec)))/60, ...
                       'recMax', max([f.rec, NaN])/60, 'socMin', f.socMin, 'hMode', f.tMode/3600, 'nCut', nnz(f.cut));
            % e: summary of the case, recovery times in minutes
            if o == 1
                R.slew.first(end + 1) = e;
            else
                R.slew.worst(end + 1) = e;
            end
        end
    end
end

% ---------------------------------------------------------------- density and decay
for i = 1:3
    [rMin, iMin] = min(rhoG(i, :)); % rMin: lowest density around the equator [kg/m^3]
    [rMax, iMax] = max(rhoG(i, :)); % rMax: highest density around the equator [kg/m^3]
    D = struct(); % D: density numbers of case i
    D.level = caseName{i};
    D.F107 = F107(i);
    D.Ap = Ap(i);
    D.rho0 = dec{i, end}.rho(1); % orbit mean at sim start [kg/m^3]
    D.min = rMin;
    D.lstMin = lstG(iMin);
    D.max = rMax;
    D.lstMax = lstG(iMax);
    D.ratio = rMax/rMin;
    D.r186 = interp1(lstG, rhoG(i, :), 18)/interp1(lstG, rhoG(i, :), 6); % 18:00 over 06:00
    R.density(i) = D;
end
R.density_h = sv{1}.h0;

nB = numel(Bg); % nB: number of ballistic coefficients
hy = NaN(3, nB, ceil(yrs)); % hy: altitude after each year, NaN once below the floor [km]
t600 = NaN(3, nB); % t600: time below hLine, NaN if never [s]
tFl = NaN(3, nB); % tFl: time below the floor, NaN if never [s]
for i = 1:3
    for j = 1:nB
        d = dec{i, j}; % d: one decay run
        hy(i, j, :) = interp1(d.t, d.h, 365.25*86400*(1:ceil(yrs)));
        k6 = find(d.h < hLine, 1); % k6: first step below hLine
        if ~isempty(k6) && k6 > 1
            t600(i, j) = d.t(k6 - 1) + (hLine - d.h(k6 - 1))/(d.h(k6) - d.h(k6 - 1))*(d.t(k6) - d.t(k6 - 1));
        elseif ~isempty(k6)
            t600(i, j) = 0;
        end
        tFl(i, j) = d.tFloor;
    end
end
h5 = hy(:, :, end); % h5: 3 x nB altitude at the end of the mission, NaN after the floor [km]
R.decay.B = Bg;
R.decay.hLine = hLine;
R.decay.hFloor = base.decay.hFloor;
R.decay.h5 = h5;
R.decay.t600 = t600/86400/365.25; % [yr]
R.decay.tFloor = tFl/86400/365.25; % [yr]

% The spacecraft of the Inputs sheet, run nB + 1:
R.decay.sheet.B = Bs;
R.decay.sheet.mass = base.drag.mass; % [kg]
R.decay.sheet.Cd = base.drag.Cd;
R.decay.sheet.area = base.drag.area; % [m^2]
for i = 1:3
    d = dec{i, nB + 1}; % d: decay run of the sheet spacecraft in case i
    R.decay.sheet.h5(i) = interp1(d.t, d.h, 365.25*86400*yrs); % altitude after the mission, NaN after the floor [km]
    k6 = find(d.h < hLine, 1); % k6: first step below hLine
    if isempty(k6)
        R.decay.sheet.t600(i) = NaN; % time below hLine, NaN if never [yr]
    elseif k6 == 1
        R.decay.sheet.t600(i) = 0;
    else
        R.decay.sheet.t600(i) = (d.t(k6 - 1) + (hLine - d.h(k6 - 1))/(d.h(k6) - d.h(k6 - 1))*(d.t(k6) - d.t(k6 - 1)))/86400/365.25;
    end
    R.decay.sheet.tFloor(i) = d.tFloor/86400/365.25; % time below the floor, NaN if never [yr]
end

% Ballistic coefficient that stays above hLine for the mission, linear in log B between the bracketing runs:
R.decay.Bstar = NaN(1, 3); % [kg/m^2]
R.decay.BstarKind = cell(1, 3); % 'below' when even the smallest B stays above, 'above' when the largest falls, 'at' otherwise
for i = 1:3
    h5i = h5(i, :); % h5i: end altitude of case i, the floor counted as zero
    h5i(isnan(h5i)) = 0;
    if h5i(1) >= hLine
        R.decay.BstarKind{i} = 'below';
        R.decay.Bstar(i) = Bg(1);
    elseif h5i(end) < hLine
        R.decay.BstarKind{i} = 'above';
        R.decay.Bstar(i) = Bg(end);
    else
        R.decay.BstarKind{i} = 'at';
        k = find(h5i >= hLine, 1); % k: first run above hLine
        R.decay.Bstar(i) = exp(log(Bg(k - 1)) + (hLine - h5i(k - 1))*(log(Bg(k)) - log(Bg(k - 1)))/(h5i(k) - h5i(k - 1)));
    end
end

% Delta-v to hold the start altitude over the mission and the starting decay rate, both at B = 1 kg/m^2 from the B = 1e6 run:
for i = 1:3
    d = dec{i, end};
    a = base.Re + d.h; % a: 1 x N semi-major axis [km]
    rate1 = abs(d.dadt)*Ball(end); % rate1: 1 x N decay rate at B = 1 kg/m^2 [km/s]
    R.decay.dvK(i) = trapz(d.t, rate1.*sqrt(base.mu./a)./(2*a))*1e3; % [m/s at B = 1 kg/m^2], dv = integral |da/dt| v/(2a) dt
    R.decay.dadtB(i) = rate1(1)*86400e3; % [m/day at B = 1 kg/m^2]
end

% Planes 2 onward against plane 1 at Bp, end altitude where both stay above the floor:
dP = []; % dP: end altitude differences [km]
for i = 1:3
    for k = 2:nPl
        for j = 1:numel(Bp)
            d = decP{i, k, j};
            d1 = dec{i, find(Bg == Bp(j), 1)}; % d1: plane 1 at the same B
            hk = interp1(d.t, d.h, 365.25*86400*yrs); % hk: end altitude of plane k [km]
            h1 = interp1(d1.t, d1.h, 365.25*86400*yrs); % h1: end altitude of plane 1 [km]
            if ~isnan(hk) && ~isnan(h1)
                dP(end + 1) = abs(hk - h1);
            end
        end
    end
end
R.decay.planeDiff = max([dP, NaN]);
R.decay.Bp = Bp;
R.decay.mass = base.drag.mass;

% ---------------------------------------------------------------- solar radiation pressure
R.srp.on = base.use.srp;
R.srp.Cr = base.srp.Cr;
R.srp.area = base.srp.area; % area facing the Sun [m^2]
R.srp.mass = base.drag.mass; % [kg]

% Size in full sunlight at 1 AU, Cr A S0/(m c):
R.srp.a1AU = base.srp.Cr*base.srp.area*base.S0/(base.drag.mass*base.c); % [m/s^2]

% ---------------------------------------------------------------- run details
R.version = '1.1.12';
R.run = datestr(now, 'yyyy-mm-dd HH:MM');
R.matlab = version;
R.figures = struct('name', {}, 'w', {}, 'h', {});

% ---------------------------------------------------------------- text summary
txtFile = fullfile(outDir, 'mission.txt'); % txtFile: printed summary
if isfile(txtFile)
    delete(txtFile);
end
diary(txtFile);
fprintf('MISSION, ODIN-Model %s, run %s, MATLAB %s\n', R.version, R.run, R.matlab);
fprintf('\norbit: %g km osculating at sim start, %.3f km orbit mean, i %g deg, nodal period %.3f min\n', base.h, R.orbit.h0, base.inc, TN/60);
fprintf('node rate %.6f deg/day against %.6f deg/day for the mean Sun: LTAN drifts %.1f min a year; sun-synchronous at this altitude is %.3f deg\n', ...
        R.orbit.dW, wSun, R.orbit.ltanDrift, R.orbit.incSSO);
fprintf('solar radiation pressure %s: Cr %g, %g m^2 facing the Sun, %g kg, %.3e m/s^2 in full sunlight at 1 AU\n', onOff(R.srp.on), R.srp.Cr, R.srp.area, R.srp.mass, R.srp.a1AU);
fprintf('\nEarth shadow over %d days from %s, drag off so the altitude holds, penumbra counted:\n', nDay, R.orbit.epoch);
for k = 1:nPl
    P = R.planes(k);
    fprintf('plane %d, LTAN %s at sim start, Sun %.1f to %.1f deg from the orbit plane: longest %.2f min on %s, beta %.2f deg, LTAN %s; no shadow on %d days, first shadow on day %d\n', ...
            k, hhmm(P.ltan0), P.betaMin, P.betaMax, P.longest, P.date, P.beta, datestr(P.ltan/24, 'HH:MM'), P.shadowFree, P.firstShadowDay);
    fprintf('   longest per year [min] %s, share of days with shadow %s\n', mat2str(round(100*P.yearLongest)/100), mat2str(round(1000*P.yearShare)/1000));
end
fprintf('worst orbit: plane %d on %s, %.2f min; first week %s to %s, worst week %s to %s\n', kW, R.worst.date, R.worst.dur, ...
        R.first.weekStart, R.first.weekEnd, R.worst.weekStart, R.worst.weekEnd);

fprintf('\nloads [W]: Nominal %.2f, Raven %.2f (%.2f during array release), Safe %.2f; components less mode loads %s W\n', pL(1), pL(2), pL(2) + ld.release, pL(3), mat2str(R.loads.sumCheck));
fprintf('currents at %.1f V: %s A, Nominal %.2f A; limits %g A per channel, %g A per pack, battery %.1f W\n', base.batt.vMin, mat2str(round(100*R.loads.iComponent)/100), ...
        R.loads.iNominal, R.limits.iChannel, R.limits.iPack, R.limits.pBatt);
fprintf('telemetry frame %d B, %.0f bit/s and %.1f MB a day at %g Hz\n', R.frame.total, R.frame.bps, R.frame.MBday, fHz);

fprintf('\nbattery: %g Wh rated, full charge %.2f Wh (socTop %.4f); through the %.2f min shadow every orbit at %.2f W for %g years:\n', base.batt.energy, eFull0, ...
        base.batt.socTop, tEcl/60, pL(1), yrs);
fprintf('   energy per shadow %.3f Wh, kR %.4f, %.0f equivalent full cycles; capacity %.4f, resistance %.3f times, discharge efficiency %.4f at the end\n', ...
        b.eDraw, b.kR, b.efc(end), b.qEOL, b.rEOL, b.etaDisEOL);
fprintf('   lowest charge, share of full charge: %.4f at the start, %.4f at the end; per year %s\n', b.socMin(1), b.socMin(end), mat2str(round(1e4*R.batt.yearSoc)/1e4));
fprintf('   lowest charge per mode at the end: %s, at the start: %s\n', mat2str(round(1e4*R.batt.socModeEOL)/1e4), mat2str(round(1e4*R.batt.socModeBOL)/1e4));

fprintf('\narray at the bus over one nodal period, %g m^2, end-of-life factor %.4f, panel from %g deg C at the start of each week:\n', base.power.panelArea, p.agingRetained, p.initialTemperature_C);
for c = 1:4
    fprintf('   %-24s peak %6.2f W, mean %6.2f W, %6.2f Wh, panel %6.1f to %6.1f deg C, %5.1f%% more with the temperature coefficient held in %g to %g deg C\n', ...
            A(c).name, A(c).peak, A(c).mean, A(c).eGen, A(c).Tmin, A(c).Tmax, 100*(A(c).eClamp/A(c).eGen - 1), p.powerTempCoeffRange_C);
end

fprintf('\nenergy margin, end-of-life battery, panels on the Sun, charge %.2f and discharge %.4f efficiency:\n', base.batt.etaChg, b.etaDisEOL);
for m = 1:3
    fprintf('   %-8s %6.2f W: worst orbit load %.3f Wh, margin %+.1f%%, battery %+.3f Wh per orbit, spare %.3f W; first week margin %+.1f%%\n', ...
            modeName{m}, pL(m), mW.eLoad(m), 100*mW.M(m), mW.net(m), mW.spare(m), 100*mS.M(m));
end
fprintf('   worst week in Nominal on the Sun: lowest charge %.4f, hours in Nominal, Raven, Safe %s\n', R.week.socMin, mat2str(round(10*R.week.hMode)/10));
fprintf('   not in the budget: two pack heaters %.2f W each, EPS electronics %.3f W\n', R.limits.heater, R.limits.epsIdle);

fprintf('\nobservation windows, each started once the battery is full again in Nominal, end-of-life battery:\n');
fprintf('%-6s %7s %7s | %8s %9s %8s | %8s %9s %8s %7s\n', 'tilt', 'window', '', 'first', 'recovery', 'lowest', 'worst', 'recovery', 'lowest', 'Raven');
for e = 1:numel(R.slew.first)
    a1 = R.slew.first(e); % a1: first week case
    a2 = R.slew.worst(e); % a2: worst week case
    fprintf('%6.1f %7d %7s | %8.2f %9.1f %7.1f%% | %8.2f %9.1f %7.1f%% %6.2f h\n', a1.tilt, a1.D, 'min', a1.perDay, a1.recMed, 100*a1.socMin, ...
            a2.perDay, a2.recMed, 100*a2.socMin, a2.hMode(2));
end

fprintf('\nNRLMSISE-00 around the equator at %.1f km on %s:\n', sv{1}.h0, R.orbit.epoch);
for i = 1:3
    D = R.density(i);
    fprintf('   %-8s F10.7 %3g Ap %2g: orbit mean %.3e, lowest %.3e at %s, highest %.3e at %s, ratio %.2f, 18:00 over 06:00 %.2f\n', D.level, D.F107, D.Ap, ...
            D.rho0, D.min, datestr(D.lstMin/24, 'HH:MM'), D.max, datestr(D.lstMax/24, 'HH:MM'), D.ratio, D.r186);
end

fprintf('\ndecay of plane 1 from %.1f km, altitude after %g years or the time below %g km:\n', sv{1}.h0, yrs, hLine);
for i = 1:3
    fprintf('   %-8s', caseName{i});
    for j = 1:nB
        if isnan(t600(i, j))
            fprintf(' %g: %.1f km;', Bg(j), h5(i, j));
        else
            fprintf(' %g: %.2f yr;', Bg(j), t600(i, j)/86400/365.25);
        end
    end
    fprintf('\n   %-8s above %g km for %g years: B %s %.2f kg/m^2; dv to hold the start %.1f/B m/s, starting decay %.2f/B m/day\n', '', hLine, yrs, ...
            R.decay.BstarKind{i}, R.decay.Bstar(i), R.decay.dvK(i), R.decay.dadtB(i));
end
fprintf('   sheet spacecraft, %g kg, Cd %g, %g m^2, B = %.2f kg/m^2:', base.drag.mass, base.drag.Cd, base.drag.area, Bs);
for i = 1:3
    if isnan(R.decay.sheet.t600(i))
        fprintf(' %s %.1f km after %g years;', caseName{i}, R.decay.sheet.h5(i), yrs);
    else
        fprintf(' %s below %g km after %.2f yr;', caseName{i}, hLine, R.decay.sheet.t600(i));
    end
end
fprintf('\n');
if nPl == 2
    fprintf('   plane 2 at B = %s kg/m^2 ends within %.2f km of plane 1\n', mat2str(Bp), R.decay.planeDiff);
else
    fprintf('   planes 2 to %d at B = %s kg/m^2 end within %.2f km of plane 1\n', nPl, mat2str(Bp), R.decay.planeDiff);
end

% ---------------------------------------------------------------- figures for the slides
fig = cell(0, 3); % fig: figure name, width and height [in]

% Longest shadow of each day over the mission, two sizes:
sz = {'shadowMission', 7.2, 4.5, 'longest shadow of each day [min]'; 'shadowAssumptions', 5.0, 3.1, 'longest shadow [min]'}; % sz: name, width [in], height [in], y label
yrD = (0:nDay - 1)/365.25; % yrD: 1 x nDay day starts [yr]
for s = 1:size(sz, 1)
    [f, ax] = slideFigure(sz{s, 1}, sz{s, 2}, sz{s, 3}, fs);
    for k = 1:nPl
        plot(ax, yrD, sv{k}.day/60, 'Color', col(k, :), 'LineWidth', 1.2);
    end
    yTop = 10*ceil(max(cellfun(@(v) v.dMax, sv))/600 + 0.5); % yTop: top of the axis, the longest shadow rounded up to 10 min with room for the legend [min]
    set(ax, 'XLim', [0 yrs], 'XTick', 0:yrs, 'YLim', [0 yTop], 'YTick', 0:10:yTop);
    xlabel(ax, sprintf('years from %s', R.orbit.epoch));
    ylabel(ax, sz{s, 4});
    legend(ax, arrayfun(@(k) sprintf('plane %d, LTAN %s', k, hhmm(R.planes(k).ltan0)), 1:nPl, 'UniformOutput', false), ...
           'Location', 'northeast', 'FontSize', fs);
    exportgraphics(f, fullfile(figDir, [sz{s, 1} '.png']), 'Resolution', 300);
    fig(end + 1, :) = sz(s, 1:3);
end

% Array power over one orbit:
[f, ax] = slideFigure('arrayOrbit', 6.6, 4.5, fs);
plot(ax, A(1).t/60, A(1).p, '--', 'Color', grey, 'LineWidth', 1.2);
plot(ax, A(3).t/60, A(3).p, 'Color', col(1, :), 'LineWidth', 1.5);
plot(ax, A(4).t/60, A(4).p, 'Color', col(2, :), 'LineWidth', 1.5);
set(ax, 'XLim', [0 TN/60], 'YLim', [0 10*ceil(max([A.peak])/10)]);
xlabel(ax, 'minutes into the orbit');
ylabel(ax, 'array power at the bus [W]');
legend(ax, {'first week, on the Sun', 'worst orbit, on the Sun', sprintf('worst orbit, %g^{\\circ} from the Sun', base.slew.tilt)}, 'Location', 'southoutside', 'FontSize', fs);
exportgraphics(f, fullfile(figDir, 'arrayOrbit.png'), 'Resolution', 300);
fig(end + 1, :) = {'arrayOrbit', 6.6, 4.5};

% Battery over the mission:
[f, ax] = slideFigure('batteryMission', 5.6, 3.6, fs);
plot(ax, b.t/365.25, 100*b.socMin, 'Color', col(1, :), 'LineWidth', 1.6);
plot(ax, b.t/365.25, 100*b.q, '--', 'Color', col(2, :), 'LineWidth', 1.4);
plot(ax, [0 yrs], 100*base.batt.floor*[1 1], 'k:', 'LineWidth', 1.0);
set(ax, 'XLim', [0 yrs], 'XTick', 0:yrs, 'YLim', [0 100], 'YTick', 0:20:100);
xlabel(ax, 'years in orbit');
ylabel(ax, 'share [%]');
legend(ax, {'lowest charge in the worst shadow', 'capacity left', sprintf('%g%% floor', 100*base.batt.floor)}, 'Location', 'southwest', 'FontSize', fs);
exportgraphics(f, fullfile(figDir, 'batteryMission.png'), 'Resolution', 300);
fig(end + 1, :) = {'batteryMission', 5.6, 3.6};

% Charge over the first two days of the worst week, Sun-pointed and with the longest windows at cfg.slew.tilt:
[f, ax] = slideFigure('windowsCharge', 11.9, 2.6, fs);
tw = (W.t - W.t(1))/86400; % tw: 1 x N days into the worst week
q = tw <= 2; % q: first two days
plot(ax, tw(q), 100*pmW{1}.soc(q), '--', 'Color', grey, 'LineWidth', 1.1);
plot(ax, tw(q), 100*frW(1, end).soc(q), 'Color', col(1, :), 'LineWidth', 1.3);
plot(ax, [0 2], 100*base.modes.socRaven*[1 1], 'k:', 'LineWidth', 1.0);
set(ax, 'XLim', [0 2], 'XTick', 0:0.5:2, 'YLim', [50 100], 'YTick', 50:10:100);
xlabel(ax, sprintf('days into the worst week, from %s', R.worst.weekStart));
ylabel(ax, 'charge [%]');
legend(ax, {'on the Sun', sprintf('%d min windows at %g^{\\circ}', slewMin(end), tilts(1)), sprintf('Raven below %g%%', 100*base.modes.socRaven)}, 'Location', 'eastoutside', 'FontSize', fs);
exportgraphics(f, fullfile(figDir, 'windowsCharge.png'), 'Resolution', 300);
fig(end + 1, :) = {'windowsCharge', 11.9, 2.6};

% Density by local solar time:
[f, ax] = slideFigure('densityLocalTime', 6.8, 3.3, fs);
for i = 1:3
    plot(ax, lstG, rhoG(i, :), 'Color', caseCol(i, :), 'LineWidth', 1.4);
end
set(ax, 'YScale', 'log', 'XLim', [0 24], 'XTick', 0:3:24);
xlabel(ax, 'local solar time at the equator [h]');
ylabel(ax, 'density [kg/m^3]');
legend(ax, strcat('ECSS', {' '}, caseName), 'Location', 'eastoutside', 'FontSize', fs);
exportgraphics(f, fullfile(figDir, 'densityLocalTime.png'), 'Resolution', 300);
fig(end + 1, :) = {'densityLocalTime', 6.8, 3.3};

% Altitude after the mission against the ballistic coefficient:
[f, ax] = slideFigure('decayB', 6.6, 4.5, fs);
for i = 1:3
    hp = h5(i, :); % hp: end altitude, runs that fell plotted at the floor [km]
    hp(isnan(hp)) = base.decay.hFloor;
    plot(ax, Bg, hp, '-o', 'Color', caseCol(i, :), 'MarkerFaceColor', caseCol(i, :), 'MarkerSize', 4, 'LineWidth', 1.3);
end
plot(ax, Bg([1 end]), hLine*[1 1], 'k:', 'LineWidth', 1.0);
yLo = base.decay.hFloor; % yLo: bottom of the axis [km]
yHi = 50*ceil((sv{1}.h0 + 10)/50); % yHi: top of the axis [km]
plot(ax, Bs*[1 1], [yLo yHi], '--', 'Color', grey, 'LineWidth', 1.2);
set(ax, 'XScale', 'log', 'XLim', Bg([1 end]), 'XTick', [10 20 50 100 200 500 1000], 'XTickLabel', {'10', '20', '50', '100', '200', '500', '1000'}, ...
    'YLim', [yLo yHi]);
xlabel(ax, 'B = m/(C_d A) [kg/m^2]');
ylabel(ax, sprintf('orbit-mean altitude after %g years [km]', yrs));
legend(ax, [strcat('ECSS', {' '}, caseName), {sprintf('%g km', hLine), sprintf('sheet, %.1f kg/m^2', Bs)}], 'Location', 'southeast', 'FontSize', fs);
exportgraphics(f, fullfile(figDir, 'decayB.png'), 'Resolution', 300);
fig(end + 1, :) = {'decayB', 6.6, 4.5};

R.figures = cell2struct(fig, {'name', 'w', 'h'}, 2).';
fprintf('\nfigures in %s: %s\n', figDir, strjoin(fig(:, 1).', ', '));

fid = fopen(fullfile(outDir, 'mission.json'), 'w'); % fid: numbers file
fprintf(fid, '%s', jsonencode(R, 'PrettyPrint', true));
fclose(fid);
fprintf('numbers in %s and %s\n', txtFile, fullfile(outDir, 'mission.json'));
diary off;


function s = onOff(tf)
% onOff: on or off for a logical switch.
% tf: switch value
% s: text for the log

if tf
    s = 'on';
else
    s = 'off';
end

end


function s = hhmm(h)
% hhmm: clock time to the nearest minute.
% h: time of day [h]
% s: text HH:MM

m = mod(round(60*h), 1440); % m: minutes after midnight
s = sprintf('%02d:%02d', floor(m/60), mod(m, 60));

end
