function dragCheck(nDay)
% dragCheck: validation of the J77 drag model, from the density code up to the decay of case A against the high-fidelity run.
% nDay: length of the case A propagations, 365 when omitted [days]

if nargin < 1
    nDay = 365;
end

here = fileparts(mfilename('fullpath')); % here: folder holding this file
root = fullfile(here, '..'); % root: repository root folder
addpath(root);
odinPath;

outDir = fullfile(here, 'results'); % outDir: folder for logs and results
if ~isfolder(outDir)
    mkdir(outDir);
end
logFile = fullfile(outDir, sprintf('dragCheck_%dd.txt', nDay)); % logFile: text log of this run
if isfile(logFile)
    delete(logFile);
end
diary(logFile);

ok = true(1, 6); % ok: pass flag of checks 1 to 6
base = config(); % base: parameter struct from config.m

Fbar = [70 100 150 200 250 300]; % Fbar: 81-day mean 10.7 cm flux of J77 Table 5 [sfu]
Ttab = [720.9 860.5 1057.2 1227.4 1380.7 1522.2]; % Ttab: J77 Table 5 exospheric temperature [K]
Tinf = zeros(size(Fbar)); % Tinf: exospheric temperature returned by dragTable.m [K]
for k = 1:numel(Fbar)
    cfg = base; % cfg: parameter struct of this check
    cfg.drag.F107 = Fbar(k);
    [~, Tinf(k)] = dragTable(cfg);
end
ok(1) = max(abs(Tinf - Ttab)) <= 0.05;
fprintf('1. exospheric temperature against J77 Table 5\n');
fprintf('   F10.7 [sfu]   Table 5 [K]   dragTable [K]\n');
fprintf('   %11.0f   %11.1f   %13.1f\n', [Fbar; Ttab; Tinf]);
fprintf('   largest difference %.3f K, limit 0.05 K: %s\n\n', max(abs(Tinf - Ttab)), verdict(ok(1)));

ref = load(fullfile(root, 'external', 'j77', 't1000.out'), '-ascii'); % ref: t1000.out rows of height [km], T [K], log10 number densities of N2, O2, O, Ar, He, H, total [1/m^3], molecular weight [g/mol]
[T, CN2, CO2, CO, CAr, CHe, CH, CM, WM] = j77sri(2500, 1000, 0:2500);
% T: temperature at 0, 1, ..., 2500 km [K]
% CN2, CO2, CO, CAr, CHe, CH, CM: number densities of N2, O2, O, Ar, He, H and their total [1/cm^3]
% WM: mean molecular weight [g/mol]
z = ref(:,1)'; % z: reference heights [km]
C = [CN2; CO2; CO; CAr; CHe; CH; CM]; % C: 7 x 2501 number densities [1/cm^3]
C = C(:, z+1);
L = -9.9 * ones(size(C)); % L: log10 number densities in the reference file's convention [log10(1/m^3)]
L(C > 1.26e-16) = log10(C(C > 1.26e-16)) + 6;
dT = max(abs(T(z+1) - ref(:,2)')); % dT: largest temperature difference [K]
dL = max(max(abs(L - ref(:,3:9)'))); % dL: largest log10 number density difference
dW = max(abs(WM(z+1) - ref(:,10)')); % dW: largest molecular weight difference [g/mol]
ok(2) = dT <= 0.006 && dL <= 6e-5 && dW <= 6e-4;
fprintf('2. j77sri.m against t1000.out, %d rows: largest differences %.4f K, %.1e in log10 n, %.4f g/mol\n', numel(z), dT, dL, dW);
fprintf('   limits 0.006 K, 6e-5, 0.0006 g/mol, just above half the last printed digit: %s\n\n', verdict(ok(2)));

cfg = base;
cfg.drag.F107 = 150;
logRho = dragTable(cfg); % logRho: J77 density interpolant at F10.7 = 150
hh = [400 500 550 600 650 700 800 900 950 1000]; % hh: comparison heights [km]
rhoJ = exp(logRho(hh)); % rhoJ: J77 density [kg/m^3]
t0 = datetime(base.epoch); % t0: UTC epoch
[LAT, LST] = ndgrid(-87.5:5:87.5, 0.5:1:23.5);
% LAT: latitude grid [deg]
% LST: local solar time grid [h]
m = numel(LAT); % m: number of grid points
lon = mod((LST(:) - hour(t0))*15 + 180, 360) - 180; % lon: longitude where the local solar time is LST at the epoch [deg]
w = cosd(LAT(:)); % w: area weight of each grid point
rhoM = zeros(size(hh)); % rhoM: NRLMSISE-00 global mean density, F10.7 = 150 and Ap = 4 [kg/m^3]
for j = 1:numel(hh)
    [~, d] = atmosnrlmsise00(hh(j)*1e3*ones(m, 1), LAT(:), lon, year(t0)*ones(m, 1), day(t0, 'dayofyear')*ones(m, 1), ...
                             zeros(m, 1), LST(:), 150*ones(m, 1), 150*ones(m, 1), 4*ones(m, 7));
    % d: m x 9 densities, column 6 the total mass density [kg/m^3]
    rhoM(j) = sum(w .* d(:,6)) / sum(w);
end
ok(3) = all(abs(rhoJ./rhoM - 1) <= 0.15);
fprintf('3. J77 at F10.7 = 150 against the NRLMSISE-00 global mean at F10.7 = 150, Ap = 4 on %s\n', string(t0, 'dd MMM yyyy'));
fprintf('   height [km]   J77 [kg/m^3]   NRLMSISE-00 [kg/m^3]   ratio\n');
fprintf('   %11.0f   %12.3e   %20.3e   %5.2f\n', [hh; rhoJ; rhoM; rhoJ./rhoM]);
fprintf('   all ratios within 15%%, the scatter of Jacchia models against drag data: %s\n\n', verdict(ok(3)));

cfg.use.drag = true;
cfg.drag.logRho = logRho;
cfg.drag.mass = 30;
cfg.drag.area = 1;
cfg.drag.Cd = 1;
cfg0 = cfg; % cfg0: the same settings with drag off
cfg0.use.drag = false;
[r, v] = oe2rv(cfg.Re + 650, 0, 97.95, 40, 0, 30, cfg.mu);
% r: test position in ECI [km]
% v: test velocity in ECI [km/s]
aD = eom(0, [r; v], cfg) - eom(0, [r; v], cfg0); % aD: drag part of the state derivative [km/s; km/s^2]

% Height above the ellipsoid, first order in flattening:
h = norm(r) - cfg.Re * (1 - cfg.f * (r(3)/norm(r))^2); % h: height above the ellipsoid [km]

% Velocity relative to an atmosphere rotating with the Earth, v - wE x r, in metres:
vRel = 1e3 * (v - cfg.wE * [-r(2); r(1); 0]); % vRel: velocity relative to the atmosphere [m/s]

% Drag acceleration in SI units:
aSI = -0.5 * exp(logRho(h)) / (cfg.drag.mass/(cfg.drag.Cd*cfg.drag.area)) * norm(vRel) * vRel; % aSI: drag acceleration [m/s^2]

eA = norm(1e3*aD(4:6) - aSI) / norm(aSI); % eA: relative difference of the two accelerations
ok(4) = eA < 1e-6 && all(aD(1:3) == 0);
fprintf('4. drag acceleration from eom.m against SI units at %.1f km: relative difference %.1e, limit 1e-6: %s\n\n', h, eA, verdict(ok(4)));

r1 = tradeLighting(1/24, 1); % r1: case A parameters from tradeLighting.m, one hour run
cfgA = r1.cfg; % cfgA: case A parameter struct, spacecraft set to m/(Cd A) = 30 kg/m^2 as in the high-fidelity run
cfgA.use.eclipse = false;
cfgA.drag.F107 = 150;
cfgA.drag.mass = 30;
cfgA.drag.area = 1;
cfgA.drag.Cd = 1;
nOrb = nDay*86400/r1.T; % nOrb: run length in orbital periods
w1 = round(r1.T/cfgA.ode.dtOut); % w1: one orbital period [samples]

tic;
out0 = propagate(cfgA, nOrb); % out0: case A without drag
tRun0 = toc; % tRun0: wall time without drag [s]

cfgA.use.drag = true;
tic;
out1 = propagate(cfgA, nOrb); % out1: case A with J77 drag
tRun1 = toc; % tRun1: wall time with drag [s]

alt0 = movmean(vecnorm(out0.r) - cfgA.Re, w1); % alt0: orbit averaged altitude without drag [km]
alt1 = movmean(vecnorm(out1.r) - cfgA.Re, w1); % alt1: orbit averaged altitude with drag [km]
k1 = numel(alt1) - w1; % k1: sample one orbit before the end
dModel = (alt1(k1) - alt1(w1)) - (alt0(k1) - alt0(w1)); % dModel: drag part of the change in orbit averaged altitude [km]

hfFile = fullfile(outDir, 'numericalCrossCheck_365d.mat'); % hfFile: saved high-fidelity results
dJ2 = NaN; % dJ2: largest difference of the drag-free run from the earlier J2 model run [km]
dHF = NaN; % dHF: drag part of the change in the high-fidelity run, NRLMSISE-00 at F10.7 = 150 and 30 kg/m^2 [km]
if isfile(hfFile)
    hf = load(hfFile, 'task', 'alt', 'altJ2'); % hf: high-fidelity configurations and orbit averaged altitudes
    kk = w1+1 : numel(alt0)-w1; % kk: samples with a full averaging window in both runs
    dJ2 = max(abs(alt0(kk) - hf.altJ2(kk, 1)'));
    jNo = find(~[hf.task.drag], 1); % jNo: high-fidelity run without drag
    jD = find([hf.task.drag] & strcmp({hf.task.sw}, 'constant') & [hf.task.B] == 30, 1); % jD: high-fidelity run with drag at F10.7 = 150 and 30 kg/m^2
    dHF = (hf.alt{jD}(k1, 1) - hf.alt{jD}(w1, 1)) - (hf.alt{jNo}(k1, 1) - hf.alt{jNo}(w1, 1));
end
ok(5) = dJ2 < 1e-6;

cfgA.drag.logRho = logRho;
a = cfgA.Re + alt0(w1); % a: semi-major axis of the orbit averaged run, started at the model's orbit averaged altitude [km]
span = (k1 - w1) * cfgA.ode.dtOut; % span: time between the two compared samples [s]
nStep = ceil(span/86400); % nStep: number of integration steps
for k = 1:nStep
    % Midpoint step of the orbit averaged decay rate:
    a = a + span/nStep * dadtCircular(a + span/nStep/2 * dadtCircular(a, cfgA), cfgA);
end
dSemi = a - cfgA.Re - alt0(w1); % dSemi: change in semi-major axis from the orbit averaged rate [km]
ok(6) = abs(dModel/dSemi - 1) <= 0.02;

if isfile(hfFile)
    fprintf('5. drag-free case A against the earlier J2 model run: largest difference %.1e km, limit 1e-6 km: %s\n', dJ2, verdict(ok(5)));
else
    fprintf('5. %s not found: %s\n', hfFile, verdict(ok(5)));
end
fprintf('6. case A drag decay over %d days, F10.7 = 150, m/(Cd A) = 30 kg/m^2\n', nDay);
fprintf('   model with J77 drag (ode45)              %8.3f km   run %.0f s with drag, %.0f s without\n', dModel, tRun1, tRun0);
fprintf('   orbit averaged rate, SI units            %8.3f km   model / rate = %.3f, limit 2%%: %s\n', dSemi, dModel/dSemi, verdict(ok(6)));
if isfile(hfFile)
    fprintf('7. high-fidelity, NRLMSISE-00 and EGM2008   %8.3f km   model / high-fidelity = %.3f, comparison only\n\n', dHF, dModel/dHF);
else
    fprintf('7. high-fidelity results not found, comparison skipped\n\n');
end
fprintf('checks passed: %d of 6\n', nnz(ok));

diary off;
save(fullfile(outDir, sprintf('dragCheck_%dd.mat', nDay)), 'ok', 'Fbar', 'Ttab', 'Tinf', 'dT', 'dL', 'dW', 'hh', 'rhoJ', 'rhoM', 'eA', ...
     'dJ2', 'dModel', 'dSemi', 'dHF', 'tRun0', 'tRun1');

end


function dadt = dadtCircular(a, cfg)
% dadtCircular: orbit averaged rate of change of semi-major axis of a circular orbit under J77 drag, in SI units.
% a: semi-major axis [km]
% cfg: parameter struct with the J77 interpolant in cfg.drag.logRho
% dadt: rate of change of semi-major axis [km/s]

u = linspace(0, 2*pi, 3601); % u: argument of latitude samples [rad]
u(end) = [];
inc = deg2rad(cfg.inc); % inc: inclination [rad]
aM = 1e3 * a; % aM: semi-major axis [m]
mu = 1e9 * cfg.mu; % mu: Earth gravitational parameter [m^3/s^2]
v = sqrt(mu/aM); % v: circular speed [m/s]

r = aM * [cos(u); cos(inc)*sin(u); sin(inc)*sin(u)]; % r: 3 x 3600 positions, node on the x axis [m]
vv = v * [-sin(u); cos(inc)*cos(u); sin(inc)*cos(u)]; % vv: 3 x 3600 velocities [m/s]

% Velocity relative to the rotating atmosphere:
vRel = vv - cfg.wE * [-r(2,:); r(1,:); zeros(size(u))]; % vRel: 3 x 3600 velocities relative to the atmosphere [m/s]

h = a - cfg.Re * (1 - cfg.f * (r(3,:)/aM).^2); % h: 1 x 3600 heights above the ellipsoid [km]
rho = exp(cfg.drag.logRho(h)); % rho: 1 x 3600 J77 density [kg/m^3]
B = cfg.drag.mass / (cfg.drag.Cd * cfg.drag.area); % B: ballistic coefficient [kg/m^2]

% Drag acceleration:
aD = -0.5 * rho / B .* vecnorm(vRel) .* vRel; % aD: 3 x 3600 drag accelerations [m/s^2]

% Semi-major axis rate from the specific energy rate, da/dt = (2 a^2/mu) <v . aD>:
dadt = 1e-3 * 2*aM^2/mu * mean(sum(vv .* aD));

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
