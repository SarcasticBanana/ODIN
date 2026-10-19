function out = propagate(cfg, nOrbits) % integrates the equations of motion with ode45.
% cfg: parameter struct from config.m
% nOrbits: number of orbital periods
% out: struct holding the state history and derived quantities

if nargin < 2
    nOrbits = 3;
end

a = cfg.Re + cfg.h; % a: semi-major axis [km]
T = 2*pi*sqrt(a^3/cfg.mu); % T: orbital period [s]

[r0, v0] = oe2rv(a, cfg.ecc, cfg.inc, cfg.raan0, cfg.argp0, cfg.nu0, cfg.mu);
% r0: initial position vector in ECI [km]
% v0: initial velocity vector in ECI [km/s]

if cfg.use.drag
    cfg.drag.logRho = dragTable(cfg);
    % cfg.drag.logRho: pchip interpolant of the natural log of J77 mass density [kg/m^3] against height above the ellipsoid [km]
end

tspan = 0 : cfg.ode.dtOut : nOrbits*T;

opts = odeset('RelTol', cfg.ode.relTol, 'AbsTol', cfg.ode.absTol);
% opts: ode45 options structure

[t, X] = ode45(@(t,x) eom(t,x,cfg), tspan, [r0; v0], opts);
% t: column vector of output times [s]
% X: N x 6 state history, one row per output time [km, km/s]

out.t = t(:)'; % out.t: 1 x N output times since sim start [s]
out.r = X(:,1:3)'; % out.r: 3 x N position history in ECI [km]
out.v = X(:,4:6)'; % out.v: 3 x N velocity history in ECI [km/s]
out.T = T; % out.T: orbital period [s]
out.a = a; % out.a: semi-major axis [km]
out.cfg = cfg; % out.cfg: copy of the parameter struct used for this run

% Specific angular momentum h = r x v
out.hvec = cross(out.r, out.v);   % out.hvec: 3 x N specific angular momentum [km^2/s]

if cfg.use.eclipse
    [out.rsun, out.sunhat] = sunVector(out.t, cfg);
    % out.rsun: 3 x N Sun position in ECI [km]
    % out.sunhat: 3 x N unit vector toward the Sun

    out.illum = shadow(out.r, out.rsun, cfg);
    % out.illum: 1 x N visible fraction of the solar disk, 0 in umbra and 1 in full sunlight
end

end
