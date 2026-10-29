% contactCheck: contactTimes edges against the toolbox elevation for ODIN 1.1 and every station in the Stations sheet over a week.

here = fileparts(mfilename('fullpath')); % here: folder holding this script
addpath(fullfile(here, '..'));
odinPath;

outDir = fullfile(here, 'results'); % outDir: folder for logs
if ~isfolder(outDir)
    mkdir(outDir);
end
logFile = fullfile(outDir, 'contactCheck.txt'); % logFile: text log of this run
if isfile(logFile)
    delete(logFile);
end
diary(logFile);

nDay = 7; % nDay: run length [days]
dt = 10; % dt: output and scenario sample spacing [s]
cfg = config(); % cfg: parameter struct of ODIN 1.1
cfg.ode.dtOut = dt;
st = stations(); % st: ground stations from the Stations sheet

% Kepler's third law:
T = 2*pi*sqrt((cfg.Re + cfg.h)^3/cfg.mu); % T: orbital period [s]

out = propagate(cfg, (nDay*86400 + 1)/T); % out: state history on the dt grid
t0 = datetime(cfg.epoch); % t0: UTC epoch
sc = satelliteScenario(t0, t0 + seconds(out.t(end)), dt); % sc: scenario on the output grid
sat = odinToScenario(sc, out, 'ODIN 1.1'); % sat: ODIN 1.1 ephemeris in GCRF

nPass = 0; % nPass: stations whose edges all agree within a second
for g = 1:numel(st)
    con = contactTimes(out, cfg, st(g)); % con: model contact windows of station g
    gs = groundStation(sc, st(g).lat, st(g).lon, 'Altitude', st(g).alt, 'MinElevationAngle', st(g).mask, 'Name', st(g).name); % gs: ground station

    x = [con.aos, con.los]; % x: 1 x 2K model edge times [s]
    hi = [con.elMax, con.elMax] >= st(g).mask + 1; % hi: 1 x 2K mask of edges of windows peaking at least 1 deg above the mask
    hi = hi(~isnan(x));
    x = x(~isnan(x));

    [~, e0] = aer(gs, sat, t0 + seconds(x)); % e0: 1 x E toolbox elevation at the model edges [deg]
    [~, eA] = aer(gs, sat, t0 + seconds(x - 0.5)); % eA: 1 x E toolbox elevation half a second earlier [deg]
    [~, eB] = aer(gs, sat, t0 + seconds(x + 0.5)); % eB: 1 x E toolbox elevation half a second later [deg]

    % Toolbox crossing of the mask minus the model edge, one Newton step from the model edge:
    dx = -(e0 - st(g).mask)./(eB - eA); % dx: 1 x E toolbox minus model edge time [s]

    [~, eM] = aer(gs, sat, t0 + seconds(con.tMax)); % eM: 1 x K toolbox elevation at the model's highest point [deg]

    fprintf('\n%s, %.3f N %.3f E, %g m, %g deg mask, %d windows over %d days\n', st(g).name, st(g).lat, st(g).lon, ...
            st(g).alt, st(g).mask, numel(con.aos), nDay);
    fprintf('toolbox minus model edge time, windows peaking at least 1 deg above the mask: median %.3f s, largest %.3f s over %d edges\n', ...
            median(abs(dx(hi))), max(abs(dx(hi))), nnz(hi));
    fprintf('the same for lower windows: largest %.3f s over %d edges\n', max([0, abs(dx(~hi))]), nnz(~hi));
    fprintf('toolbox minus model highest elevation: largest %.4f deg\n', max(abs(eM - con.elMax)));
    nPass = nPass + all(abs(dx(hi)) < 1);
end

fprintf('\nchecks passed: %d of %d\n', nPass, numel(st));
diary off;
