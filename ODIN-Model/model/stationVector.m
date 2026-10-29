function [rgs, up] = stationVector(t, st, cfg)
% stationVector: ground station position and local vertical in the model frame over the propagation window.
% t: 1 x N times since sim start [s]
% st: station struct from stations.m
% cfg: parameter struct from config.m
% rgs: 3 x N station position [km]
% up: 3 x N unit normal to the ellipsoid at the station

lat = deg2rad(st.lat); % lat: geodetic latitude [rad]
h = st.alt/1e3; % h: height above the ellipsoid [km]
e2 = cfg.f*(2 - cfg.f); % e2: squared eccentricity of the ellipsoid

% Prime vertical radius of curvature of the ellipsoid:
N = cfg.Re/sqrt(1 - e2*sin(lat)^2); % N: prime vertical radius of curvature [km]

lam = earthRotation(t, cfg) + deg2rad(st.lon); % lam: 1 x N right ascension of the station meridian [rad]

% Geodetic to Cartesian position on the rotating ellipsoid:
rgs = [(N + h)*cos(lat)*cos(lam)
       (N + h)*cos(lat)*sin(lam)
       (N*(1 - e2) + h)*sin(lat)*ones(size(lam))];

% Ellipsoid normal at the geodetic latitude:
up = [cos(lat)*cos(lam)
      cos(lat)*sin(lam)
      sin(lat)*ones(size(lam))];

end
