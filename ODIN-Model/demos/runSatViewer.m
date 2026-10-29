addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
odinPath;
cfg = config(); % cfg: parameter struct of satellite 1.1 from config.m
% Kepler's third law:
T = 2*pi*sqrt((cfg.Re + cfg.h)^3/cfg.mu); % T: orbital period [s]
out = propagate(cfg, 86400/T); % out: state history over one day from propagate.m
t0 = datetime(cfg.epoch);
sc = satelliteScenario(t0, t0 + seconds(out.t(end)), 60);
sat = odinToScenario(sc, out, 'ODIN 1.1'); % sat: Satellite object of ODIN 1.1
st = stations(); % st: ground stations from the Stations sheet
for g = 1:numel(st)
    gs = groundStation(sc, st(g).lat, st(g).lon, 'Altitude', st(g).alt, 'MinElevationAngle', st(g).mask, 'Name', st(g).name); % gs: ground station
    access(sat, gs);
end
groundTrack(sat, 'LeadTime', 5400, 'TrailTime', 5400);
satelliteScenarioViewer(sc);
play(sc)