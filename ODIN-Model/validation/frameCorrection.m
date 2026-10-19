% frameCorrection: size of the mean of epoch to GCRF rotation applied by odinToScenario, per trade case over one day.

here = fileparts(mfilename('fullpath')); % here: folder holding this script
addpath(fullfile(here, '..'));
odinPath;

outDir = fullfile(here, 'results'); % outDir: folder for logs and figures
if ~isfolder(outDir)
    mkdir(outDir);
end
logFile = fullfile(outDir, 'frameCorrection.txt'); % logFile: text log of this run
if isfile(logFile)
    delete(logFile);
end
diary(logFile);

nDay = 1; % nDay: run length [days]
nCase = 6; % nCase: number of trade cases

dr = zeros(nCase, 3); % dr: minimum, mean and maximum position change from the rotation [km]
dA = zeros(nCase, 2); % dA: largest along-track change [km] and its time equivalent [s]
dC = zeros(nCase, 1); % dC: largest cross-track change [km]
dg = zeros(nCase, 3); % dg: minimum, mean and maximum ground track shift [km]
dI = zeros(nCase, 2); % dI: largest 60 s ephemeris interpolation error at 1 s sampling, and largest error at the ephemeris samples [m]
lbl = cell(nCase, 1); % lbl: case labels

for k = 1:nCase
    res = tradeLighting(nDay, k); % res: case parameters and shadow times from tradeLighting.m
    cfg = res.cfg; % cfg: parameter struct of this case
    lbl{k} = res.label;

    out = propagate(cfg, nDay*86400/res.T); % out: state history on the cfg.ode.dtOut grid
    t0 = datetime(cfg.epoch); % t0: UTC epoch
    tt = t0 + seconds(out.t(:)); % tt: N x 1 UTC sample times

    sc = satelliteScenario(t0, tt(end), cfg.ode.dtOut); % sc: scenario sampled on the output grid
    [sat, R] = odinToScenario(sc, out, 'corrected'); % sat: satellite with the rotation applied
    satU = satellite(sc, timetable(tt, 1e3*out.r.'), timetable(tt, 1e3*out.v.'), ...
                     'CoordinateFrame', 'inertial', 'Name', 'uncorrected'); % satU: satellite with the rotation omitted

    d = R*out.r - out.r; % d: 3 x N position change from the rotation [km]
    uA = out.v ./ vecnorm(out.v); % uA: 3 x N along-track unit vectors
    uC = out.hvec ./ vecnorm(out.hvec); % uC: 3 x N orbit normal unit vectors

    dr(k,:) = [min(vecnorm(d)), mean(vecnorm(d)), max(vecnorm(d))];
    [dA(k,1), m] = max(abs(dot(d, uA, 1))); % m: sample of the largest along-track change

    % Along-track distance over orbital speed:
    dA(k,2) = dA(k,1) / norm(out.v(:,m));

    dC(k) = max(abs(dot(d, uC, 1)));

    llaC = states(sat, 'CoordinateFrame', 'geographic'); % llaC: 3 x N latitude [deg], longitude [deg] and altitude [m] with the rotation
    llaU = states(satU, 'CoordinateFrame', 'geographic'); % llaU: 3 x N latitude [deg], longitude [deg] and altitude [m] without it

    % Haversine distance between the two subsatellite points on a sphere of radius cfg.Re:
    g = 2*cfg.Re*asin(sqrt(sind((llaC(1,:) - llaU(1,:))/2).^2 + ...
        cosd(llaC(1,:)).*cosd(llaU(1,:)).*sind((llaC(2,:) - llaU(2,:))/2).^2)); % g: 1 x N ground track shift [km]
    dg(k,:) = [min(g), mean(g), max(g)];

    cfg1 = cfg; % cfg1: parameter struct with a 1 s output grid
    cfg1.ode.dtOut = 1;
    out1 = propagate(cfg1, 2); % out1: state history over two orbits on a 1 s grid
    sc1 = satelliteScenario(t0, t0 + seconds(out1.t(end)), 1); % sc1: scenario sampled every second
    s60 = odinToScenario(sc1, out, 'eph60'); % s60: satellite from the 60 s ephemeris
    s1 = odinToScenario(sc1, out1, 'eph1'); % s1: satellite from the 1 s ephemeris
    p60 = states(s60); % p60: 3 x M GCRF position interpolated from the 60 s ephemeris [m]
    p1 = states(s1); % p1: 3 x M GCRF position from the 1 s ephemeris [m]
    dI(k,:) = [max(vecnorm(p60 - p1)), max(vecnorm(p1 - 1e3*R*out1.r))];
end

fprintf('\nmean equator and equinox of %s UTC to GCRF\n', string(t0, 'yyyy-MM-dd HH:mm'));

% Rotation angle from the trace, cos(angle) = (trace(R) - 1)/2:
fprintf('total rotation angle %.4f deg (%.1f arcsec)\n', acosd((trace(R) - 1)/2), 3600*acosd((trace(R) - 1)/2));

fprintf('equinox of epoch in GCRF: right ascension %.4f deg, declination %.4f deg\n', atan2d(R(2,1), R(1,1)), asind(R(3,1)));
fprintf('pole of epoch from the GCRF pole: %.4f deg\n', acosd(R(3,3)));

fprintf('\ncase  |dr| min  mean   max [km]   along max [km]   [s]   cross max [km]   ground shift min  mean   max [km]   60 s eph at 1 s [m]   at samples [m]\n');
for k = 1:nCase
    fprintf('%-4s  %8.1f %5.1f %6.1f   %14.1f %5.2f   %14.1f   %16.1f %5.1f %6.1f   %19.3f   %14.2e\n', ...
            lbl{k}(1), dr(k,:), dA(k,:), dC(k), dg(k,:), dI(k,:));
end

diary off;
save(fullfile(outDir, 'frameCorrection.mat'), 'dr', 'dA', 'dC', 'dg', 'dI', 'lbl', 'R');
