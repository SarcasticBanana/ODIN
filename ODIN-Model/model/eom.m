function dxdt = eom(t, x, cfg)
% eom: state derivative for ode45 under two-body gravity plus J2 and atmospheric drag.
% t: time elapsed [s]
% x: state vector [rx; ry; rz; vx; vy; vz] in ECI [km; km/s]
% cfg: parameter struct from config.m
% dxdt: state derivative [vx; vy; vz; ax; ay; az] [km/s; km/s^2]

r = x(1:3); % r: position vector in ECI [km]
v = x(4:6); % v: velocity vector in ECI [km/s]
rn = norm(r); % rn: orbital radius  [km]

% Two-body acceleration, Newton's law of gravitation for a point mass:
a = -cfg.mu / rn^3 * r; % a: total inertial acceleration [km/s^2]

if cfg.use.J2
    a = a + accelJ2(r, cfg);
end

if cfg.use.drag
    a = a + accelDrag(r, v, cfg);
end

dxdt = [v; a];

end


function aJ2 = accelJ2(r, cfg)
% accelJ2: perturbing acceleration from the J2 zonal harmonic.
% r: position vector in ECI [km]
% cfg: parameter struct from config.m
% aJ2: J2 perturbing acceleration in ECI [km/s^2]

x = r(1); % x: ECI x component of position [km]
y = r(2); % y: ECI y component of position [km]
z = r(3); % z: ECI z component of position [km]
rn = norm(r); % rn: orbital radius [km]

k = -1.5 * cfg.J2 * cfg.mu * cfg.Re^2 / rn^5; % k: common scalar factor [1/s^2]
zr2 = (z/rn)^2; % zr2: squared sine of geocentric latitude

% J2 acceleration, the gradient of the degree-2 zonal term of the
% geopotential U = -(mu/|r|) * [1 - J2*(Re/|r|)^2 * (3*sin(lat)^2 - 1)/2]:

aJ2 = k * [ x * (1 - 5*zr2)
            y * (1 - 5*zr2)
            z * (3 - 5*zr2) ];

end


function aD = accelDrag(r, v, cfg)
% accelDrag: perturbing acceleration from drag in a J77 atmosphere rotating with the Earth.
% r: position vector in ECI [km]
% v: velocity vector in ECI [km/s]
% cfg: parameter struct from config.m, with cfg.drag.logRho set by propagate.m
% aD: drag acceleration in ECI [km/s^2]

rn = norm(r); % rn: orbital radius [km]

% Height above the WGS84 ellipsoid, first order in flattening:
h = rn - cfg.Re * (1 - cfg.f * (r(3)/rn)^2); % h: height above the ellipsoid [km]

rho = exp(cfg.drag.logRho(h)); % rho: J77 mass density [kg/m^3]

if isnan(rho)
    error('eom:dragHeight', 'height %.1f km is outside the J77 density table', h);
end

% Velocity relative to an atmosphere rotating rigidly with the Earth, v - wE x r:
vRel = v - cfg.wE * [-r(2); r(1); 0]; % vRel: velocity relative to the atmosphere [km/s]

% Drag acceleration, 1e3 converting (kg/m^3)(m^2/kg)(km/s)^2 to km/s^2:
aD = -0.5e3 * rho * cfg.drag.Cd * cfg.drag.area / cfg.drag.mass * norm(vRel) * vRel;

end
