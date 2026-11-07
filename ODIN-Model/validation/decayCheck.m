function decayCheck(nDay)
% decayCheck: validation of the NRLMSISE-00 drag, atmosnrlmsise00 against pymsis 0.13 reference densities, then dragDecay.m against ode45 with NRLMSISE-00 at every step for the ECSS low, moderate and high cases.
% nDay: length of the ode45 comparison, 30 when omitted [days]

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
logFile = fullfile(outDir, sprintf('decayCheck_%dd.txt', nDay)); % logFile: text log of this run
if isfile(logFile)
    delete(logFile);
end
diary(logFile);

ok = true(1, 4); % ok: pass flag of checks 1 to 4
base = config(); % base: parameter struct from config.m
base.drag.model = 'NRLMSISE-00';

P = [400  -54.58    50.37  2026  134  32011.2   65   65   0  5.752708e-13
     500  -24.66   104.59  2027  171  78204.4  140  140  15  3.241496e-13
     550  -54.85    55.00  2028   36  83545.5  250  250  45  1.177067e-12
     600   71.37    48.91  2029  109  65036.1   65   65   0  1.317091e-14
     625    2.58   117.32  2030  240  29273.4  140  140  15  2.295153e-13
     650  -37.76   -98.52  2026  164  45430.6  250  250  45  1.284878e-13
     700  -11.74    58.75  2027  304  38681.4   65   65   0  8.595550e-15
     800  -22.92  -109.66  2028    5  51396.4  140  140  15  8.220675e-15
     900  -11.00   -72.00  2029  353  75567.5  250  250  45  1.269539e-13
     625   50.57    38.42  2030   77  29816.7   65   65   0  1.324578e-14
     625   75.96    22.82  2026  218  77798.8  140  140  15  1.020847e-13
     625  -30.71    70.56  2027  158  27114.1  250  250  45  9.632059e-13];
% P: reference points, height [km], geodetic latitude and longitude [deg], year, day of year, UTC time of day [s], F10.7 81-day mean and daily [sfu], daily Ap, and pymsis 0.13 NRLMSISE-00 total mass density with anomalous oxygen [kg/m^3]

m = size(P, 1); % m: number of reference points

% Local solar time, UTC hours plus longitude over 15 deg/h, as in dragDensity.m:
lst = P(:,6)/3600 + P(:,3)/15; % lst: m x 1 local solar time [h]

[~, d] = atmosnrlmsise00(1e3*P(:,1), P(:,2), P(:,3), P(:,4), P(:,5), P(:,6), lst, P(:,7), P(:,8), P(:,9)*ones(1, 7), ...
                         ones(1, 23), 'Oxygen', 'None');
% d: m x 9 densities, column 6 the total mass density with anomalous oxygen [kg/m^3]

e1 = max(abs(d(:,6)./P(:,10) - 1)); % e1: largest relative difference from the reference densities
ok(1) = e1 <= 1e-4;
fprintf('1. atmosnrlmsise00 with anomalous oxygen against pymsis 0.13 at %d points, 400 to 900 km, F10.7 65 to 250, Ap 0 to 45\n', m);
fprintf('   largest relative difference %.1e, limit 1e-4: %s\n\n', e1, verdict(ok(1)));

F107 = [base.decay.F107Low, base.decay.F107Moderate, base.decay.F107High]; % F107: 10.7 cm flux of the ECSS low, moderate and high cases [sfu]
Ap = [base.decay.ApLow, base.decay.ApModerate, base.decay.ApHigh]; % Ap: daily geomagnetic index of the same cases
name = {'low', 'moderate', 'high'}; % name: ECSS activity level of each case

cA = base; % cA: plane 1 without drag
cA.use.eclipse = false;
cA.use.drag = false;
cA.use.srp = false; % drag only, as in dragDecay.m

% Kepler's third law:
T = 2*pi*sqrt((cA.Re + cA.h)^3/cA.mu); % T: orbital period [s]

nOrb = nDay*86400/T; % nOrb: run length in orbital periods
w1 = round(T/cA.ode.dtOut); % w1: one orbital period [samples]

tic;
out0 = propagate(cA, nOrb); % out0: plane 1 without drag
tRun0 = toc; % tRun0: wall time without drag [s]

alt0 = movmean(vecnorm(out0.r) - cA.Re, w1); % alt0: orbit averaged altitude without drag [km]
k1 = numel(alt0) - w1; % k1: sample one orbit before the end
span = (k1 - w1)*cA.ode.dtOut; % span: time between the two compared samples [s]

r = out0.r(:, w1); % r: position one orbit after the start [km]
v = out0.v(:, w1); % v: velocity one orbit after the start [km/s]
hv = cross(r, v); % hv: specific angular momentum [km^2/s]
raanS = mod(atan2d(hv(1), -hv(2)), 360); % raanS: right ascension of the ascending node one orbit after the start [deg]
nd = [cosd(raanS); sind(raanS); 0]; % nd: unit vector to the ascending node

% Argument of latitude, the angle from the node to r about the orbit normal:
uS = atan2d(dot(cross(hv/norm(hv), nd), r), dot(nd, r)); % uS: argument of latitude one orbit after the start [deg]

dModel = zeros(1, 3); % dModel: drag part of the change in orbit averaged altitude from ode45 [km]
dSemi = zeros(1, 3); % dSemi: change in semi-major axis from dragDecay.m [km]
tRun1 = zeros(1, 3); % tRun1: wall time of each ode45 run with drag [s]
for i = 1:3
    c = cA; % c: plane 1 with NRLMSISE-00 drag in case i
    c.use.drag = true;
    c.drag.F107 = F107(i);
    c.drag.Ap = Ap(i);

    tic;
    out1 = propagate(c, nOrb); % out1: plane 1 with drag
    tRun1(i) = toc;

    alt1 = movmean(vecnorm(out1.r) - c.Re, w1); % alt1: orbit averaged altitude with drag [km]
    dModel(i) = (alt1(k1) - alt1(w1)) - (alt0(k1) - alt0(w1));

    cD = c; % cD: the orbit averaged run, started one orbit after the start at the drag-free orbit averaged altitude
    cD.h = alt0(w1);
    cD.raan0 = raanS;
    cD.argp0 = uS;
    cD.nu0 = 0;
    cD.epoch = datevec(datenum(c.epoch) + out0.t(w1)/86400);
    cD.decay.dt = span/86400/ceil(span/86400/c.decay.dt);
    dec = dragDecay(cD, span/86400); % dec: orbit averaged decay over the same span
    dSemi(i) = interp1(dec.t, dec.h, span) - dec.h(1);

    ok(i + 1) = abs(dModel(i)/dSemi(i) - 1) <= 0.02;
end

fprintf('2 to 4. plane 1 drag decay over %d days, NRLMSISE-00, m/(Cd A) = %.1f kg/m^2, drag-free run %.0f s\n', nDay, ...
        cA.drag.mass/(cA.drag.Cd*cA.drag.area), tRun0);
for i = 1:3
    fprintf('   ECSS %-8s ode45 %9.4f km (run %.0f s), dragDecay %9.4f km, ratio %.4f, limit 2%%: %s\n', name{i}, dModel(i), tRun1(i), ...
            dSemi(i), dModel(i)/dSemi(i), verdict(ok(i + 1)));
end
fprintf('\nchecks passed: %d of 4\n', nnz(ok));

diary off;
save(fullfile(outDir, sprintf('decayCheck_%dd.mat', nDay)), 'ok', 'e1', 'dModel', 'dSemi', 'tRun0', 'tRun1');

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
