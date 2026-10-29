function [cfgs, names] = constellation(cfg)
% constellation: parameter struct of every satellite, planes stepped in RAAN and the satellites of a plane stepped along one trajectory.
% cfg: parameter struct from config.m, its orbit elements those of satellite 1.1
% cfgs: 1 x N parameter structs, plane by plane, N = cfg.const.nPlanes*cfg.const.nPerPlane
% names: 1 x N satellite names, 'ODIN plane.slot'

nP = cfg.const.nPlanes; % nP: number of orbit planes
nS = cfg.const.nPerPlane; % nS: satellites per plane
cfgs = repmat(cfg, 1, nP*nS);
names = cell(1, nP*nS);
opts = odeset('RelTol', cfg.ode.relTol, 'AbsTol', cfg.ode.absTol); % opts: ode45 options structure

for k = 1:nP
    c1 = cfg; % c1: parameter struct of satellite k.1

    % RAAN stepped by cfg.const.dRaan per plane:
    c1.raan0 = mod(cfg.raan0 + (k - 1)*cfg.const.dRaan, 360);

    Tn = nodalPeriod(c1); % Tn: time between ascending node crossings [s]
    [r0, v0] = oe2rv(c1.Re + c1.h, c1.ecc, c1.inc, c1.raan0, c1.argp0, c1.nu0, c1.mu);
    % r0: initial position of satellite k.1 in ECI [km]
    % v0: initial velocity of satellite k.1 in ECI [km/s]

    cI = c1; % cI: parameter struct for the integration, carrying the drag table when drag is on
    if cI.use.drag
        cI.drag.logRho = dragTable(cI);
    end

    for j = 1:nS
        n = (k - 1)*nS + j; % n: satellite index, plane by plane
        cfgs(n) = c1;
        names{n} = sprintf('ODIN %d.%d', k, j);

        % Along-track step as a share of the nodal period:
        dt = (j - 1)*cfg.const.dU/360*Tn; % dt: time from satellite k.1 along its trajectory [s]

        if dt > 0
            [~, X] = ode45(@(t, x) eom(t, x, cI), [0 dt], [r0; v0], opts);
            % X: state history of satellite k.1 up to dt, last row the state of satellite k.j [km, km/s]

            [a, cfgs(n).ecc, cfgs(n).inc, cfgs(n).raan0, cfgs(n).argp0, cfgs(n).nu0] = rv2oe(X(end,1:3)', X(end,4:6)', cfg.mu);
            % a: osculating semi-major axis of satellite k.j [km]
            cfgs(n).h = a - cfg.Re;
        end
    end
end

end


function Tn = nodalPeriod(cfg)
% nodalPeriod: time between the first two ascending node crossings of the propagated orbit.
% cfg: parameter struct of one satellite
% Tn: nodal period [s]

cfg.use.eclipse = false;
cfg.ode.dtOut = 10;
out = propagate(cfg, 2.2); % out: state history over 2.2 orbital periods from propagate.m
z = out.r(3,:); % z: 1 x N ECI z component of position [km]
k = find(z(1:end-1) <= 0 & z(2:end) > 0, 2); % k: samples preceding the first two ascending node crossings

% Crossing times by linear interpolation of z between samples:
tc = out.t(k) - z(k).*(out.t(k+1) - out.t(k))./(z(k+1) - z(k)); % tc: ascending node crossing times [s]
Tn = tc(2) - tc(1);

end


function [a, ecc, incDeg, raanDeg, argpDeg, nuDeg] = rv2oe(r, v, mu)
% rv2oe: classical orbital elements from ECI position and velocity, the inverse of oe2rv.m.
% r: position vector in ECI [km]
% v: velocity vector in ECI [km/s]
% mu: gravitational parameter [km^3/s^2]
% a: semi-major axis [km]
% ecc: eccentricity
% incDeg: inclination [deg]
% raanDeg: right ascension of the ascending node [deg]
% argpDeg: argument of perigee [deg]
% nuDeg: true anomaly [deg]

h = cross(r, v); % h: specific angular momentum [km^2/s]
nd = cross([0; 0; 1], h); % nd: node vector [km^2/s]

% Vis-viva:
a = 1 / (2/norm(r) - dot(v, v)/mu);

% Eccentricity vector from the angular momentum, e = v x h/mu - r/|r|:
e = cross(v, h)/mu - r/norm(r); % e: eccentricity vector
ecc = norm(e);

% Inclination from h, RAAN from the node vector:
incDeg = acosd(h(3)/norm(h));
raanDeg = mod(atan2d(nd(2), nd(1)), 360);

% Argument of latitude, the angle from the node to r in the orbit plane:
u = atan2d(dot(cross(nd, r), h)/norm(h), dot(nd, r)); % u: argument of latitude [deg]

% Argument of perigee, the angle from the node to e in the orbit plane:
argpDeg = mod(atan2d(dot(cross(nd, e), h)/norm(h), dot(nd, e)), 360);

% True anomaly from the argument of latitude:
nuDeg = mod(u - argpDeg, 360);

end
