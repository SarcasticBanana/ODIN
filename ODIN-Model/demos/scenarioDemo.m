function v = scenarioDemo(idx, frame, nHours, basemap, tracks)
% scenarioDemo: satellite scenario viewer of the ODIN constellation propagated by propagate.m, for meetings.
% idx: satellites to show as indices into constellation.m, plane by plane, all when omitted; a single index makes the camera follow that satellite
% frame: camera frame, 'ECEF' when omitted, or 'Inertial' with a Sun direction marker
% nHours: scenario length from cfg.epoch, 24 when omitted [h]
% basemap: globe imagery, 'satellite' when omitted; 'darkwater' needs no internet
% tracks: true to draw each satellite's previous and next orbit as ground tracks in the ECEF frame, false when omitted
% v: satelliteScenarioViewer object

addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
odinPath;

persistent keep % keep: scenario of the latest call, held so its viewer outlives this function

if nargin < 1
    idx = [];
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
if nargin < 5 || isempty(tracks)
    tracks = false;
end

frame = validatestring(frame, {'ECEF', 'Inertial'});

dt = 10; % dt: output and scenario sample spacing [s]
speed = 300; % speed: playback speed multiplier
col = [0.00 0.45 0.70
       0.84 0.37 0.00
       0.00 0.62 0.45
       0.80 0.47 0.65
       0.90 0.62 0.00
       0.34 0.71 0.91]; % col: 6 x 3 plane colours, rows matching the trade study figures

base = config(); % base: parameter struct from config.m
[cfgs, names] = constellation(base);
% cfgs: parameter struct of each satellite from constellation.m
% names: name of each satellite
if isempty(idx)
    idx = 1:numel(cfgs);
end
nS = base.const.nPerPlane; % nS: satellites per plane
w = 0.6*(0:nS-1)/max(nS - 1, 1); % w: share of white in the shade of each slot

t0 = datetime(base.epoch); % t0: UTC epoch
sc = satelliteScenario(t0, t0 + hours(nHours), dt); % sc: scenario over the demo window
keep = sc;
st = stations(); % st: ground stations from the Stations sheet
gs = cell(1, numel(st)); % gs: GroundStation object of each station
for g = 1:numel(st)
    gs{g} = groundStation(sc, st(g).lat, st(g).lon, 'Altitude', st(g).alt, 'MinElevationAngle', st(g).mask, 'Name', st(g).name);
end

nSat = numel(idx); % nSat: number of satellites shown
sat = cell(1, nSat); % sat: Satellite object of each satellite shown
T = zeros(1, nSat); % T: orbital period of each satellite shown [s]
c = zeros(nSat, 3); % c: colour of each satellite shown

for j = 1:nSat
    k = ceil(idx(j)/nS); % k: plane of this satellite
    s = idx(j) - (k - 1)*nS; % s: slot of this satellite in its plane

    % Plane colour mixed with white for the slot shade:
    c(j,:) = (1 - w(s))*col(k,:) + w(s);

    cfg = cfgs(idx(j)); % cfg: parameter struct of this satellite
    cfg.ode.dtOut = dt;

    % Kepler's third law:
    T(j) = 2*pi*sqrt((cfg.Re + cfg.h)^3/cfg.mu);

    out = propagate(cfg, (nHours*3600 + 1)/T(j)); % out: state history on the dt grid, one second past the window

    sat{j} = odinToScenario(sc, out, names{idx(j)});
    sat{j}.MarkerColor = c(j,:);
    sat{j}.LabelFontColor = c(j,:);
    sat{j}.Orbit.LineColor = c(j,:);

    for g = 1:numel(gs)
        ac = access(sat{j}, gs{g}); % ac: line of sight from the station, drawn while the satellite is in view
        ac.LineColor = c(j,:);
    end
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
    if strcmp(frame, 'ECEF') && tracks
        groundTrack(sat{j}, 'LeadTime', T(j), 'TrailTime', T(j), 'LeadLineColor', c(j,:), ...
                    'TrailLineColor', 0.6*c(j,:), 'Viewer', v);
    end
end

if nSat == 1
    camtarget(v, sat{1});
end

play(v);

end
