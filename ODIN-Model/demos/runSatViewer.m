addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
odinPath;
r = tradeLighting(1/24, 5);                    % case B parameters
cfg = r.cfg;
out = propagate(cfg, 86400/r.T);               % one day at 60 s
t0 = datetime(cfg.epoch);
sc = satelliteScenario(t0, t0 + seconds(out.t(end)), 60);
sat = odinToScenario(sc, out, 'ODIN');
gs = groundStation(sc, 40.007, -105.266, 'MinElevationAngle', 10, 'Name', 'Boulder placeholder');
access(sat, gs);                               % draws a line to the station while it's in view
groundTrack(sat, 'LeadTime', 5400, 'TrailTime', 5400);
satelliteScenarioViewer(sc);
play(sc)