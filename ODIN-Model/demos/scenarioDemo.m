function v = scenarioDemo(cases, frame, nHours, basemap)
% scenarioDemo: satellite scenario viewer of the ODIN trade orbits propagated by propagate.m, for meetings.
% cases: trade case letters to show, 'ABCDEF' when omitted; a single letter makes the camera follow that satellite
% frame: camera frame, 'ECEF' with ground tracks when omitted, or 'Inertial' with a Sun direction marker
% nHours: scenario length from cfg.epoch, 24 when omitted [h]
% basemap: globe imagery, 'satellite' when omitted; 'darkwater' needs no internet
% v: satelliteScenarioViewer object

addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
odinPath;

persistent keep % keep: scenario of the latest call, held so its viewer outlives this function

if nargin < 1 || isempty(cases)
    cases = 'ABCDEF';
end
if nargin < 2 || isempty(frame)
    frame = 'ECEF';
end
if nargin < 3 || isempty(nHours)
    nHours = 24;
end
if nargin < 4 || isempty(basemap)
    basemap = 'satellite';
end

cases = upper(char(cases));
frame = validatestring(frame, {'ECEF', 'Inertial'});

dt = 10; % dt: output and scenario sample spacing [s]
speed = 300; % speed: playback speed multiplier
col = [0.00 0.45 0.70
       0.84 0.37 0.00
       0.00 0.62 0.45
       0.80 0.47 0.65
       0.90 0.62 0.00
       0.34 0.71 0.91]; % col: 6 x 3 case colors matching the trade study figures

base = config(); % base: parameter struct from config.m
t0 = datetime(base.epoch); % t0: UTC epoch
sc = satelliteScenario(t0, t0 + hours(nHours), dt); % sc: scenario over the demo window
keep = sc;
gs = groundStation(sc, 40.007, -105.266, 'MinElevationAngle', 10, 'Name', 'Boulder placeholder'); % gs: ground station

nSat = numel(cases); % nSat: number of trade cases shown
sat = cell(1, nSat); % sat: Satellite object of each case shown
T = zeros(1, nSat); % T: orbital period of each case shown [s]
k = zeros(1, nSat); % k: trade case index of each case shown

for j = 1:nSat
    k(j) = double(cases(j)) - double('A') + 1;
    r1 = tradeLighting(1/24, k(j)); % r1: case parameters from tradeLighting.m, one hour run
    cfg = r1.cfg; % cfg: parameter struct of this case
    cfg.ode.dtOut = dt;
    out = propagate(cfg, (nHours*3600 + 1)/r1.T); % out: state history on the dt grid, one second past the window
    T(j) = r1.T;

    sat{j} = odinToScenario(sc, out, ['ODIN ' r1.name]);
    sat{j}.MarkerColor = col(k(j),:);
    sat{j}.LabelFontColor = col(k(j),:);
    sat{j}.Orbit.LineColor = col(k(j),:);

    ac = access(sat{j}, gs); % ac: line of sight from the station, drawn while the satellite is in view
    ac.LineColor = col(k(j),:);
end

if strcmp(frame, 'Inertial')
    [~, sHat] = sunVector(out.t, cfg); % sHat: 3 x N unit vector toward the Sun

    mk.t = out.t; % mk.t: 1 x N sample times [s]

    % Marker three Earth radii along the Sun direction:
    mk.r = 3*cfg.Re*sHat; % mk.r: 3 x N marker position [km]

    % Marker velocity by central differences:
    mk.v = gradient(mk.r, dt); % mk.v: 3 x N marker velocity [km/s]

    mk.cfg = cfg; % mk.cfg: parameter struct carrying cfg.epoch

    sun = odinToScenario(sc, mk, 'Sun direction'); % sun: Sun direction marker
    sun.MarkerColor = [1 0.85 0];
    sun.LabelFontColor = [1 0.85 0];
    sun.Orbit.LineColor = [1 0.85 0];
end

v = satelliteScenarioViewer(sc, 'CameraReferenceFrame', frame, 'Basemap', basemap, 'PlaybackSpeedMultiplier', speed);

for j = 1:nSat
    if strcmp(frame, 'ECEF')
        groundTrack(sat{j}, 'LeadTime', T(j), 'TrailTime', T(j), 'LeadLineColor', col(k(j),:), ...
                    'TrailLineColor', 0.6*col(k(j),:), 'Viewer', v);
    end
end

if nSat == 1
    camtarget(v, sat{1});
end

play(v);

end
