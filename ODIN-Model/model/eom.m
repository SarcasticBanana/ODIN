function dxdt = eom(t, x, cfg)
% eom: state derivative for ode45 under two-body gravity plus J2, atmospheric drag and solar radiation pressure.
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
    a = a + accelDrag(t, r, v, cfg);
end

if cfg.use.srp
    a = a + accelSRP(t, r, cfg);
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


function aD = accelDrag(t, r, v, cfg)
% accelDrag: perturbing acceleration from drag in an atmosphere rotating with the Earth, density from dragDensity.m.
% t: time elapsed [s]
% r: position vector in ECI [km]
% v: velocity vector in ECI [km/s]
% cfg: parameter struct from config.m, with cfg.drag.logRho set by propagate.m for J77
% aD: drag acceleration in ECI [km/s^2]

rho = dragDensity(t, r, cfg); % rho: mass density [kg/m^3]

if isnan(rho)
    rn = norm(r); % rn: orbital radius [km]

    % Height above the WGS84 ellipsoid, first order in flattening:
    h = rn - cfg.Re * (1 - cfg.f * (r(3)/rn)^2); % h: height above the ellipsoid [km]

    error('eom:dragHeight', 'height %.1f km is outside the J77 density table', h);
end

% Velocity relative to an atmosphere rotating rigidly with the Earth, v - wE x r:
vRel = v - cfg.wE * [-r(2); r(1); 0]; % vRel: velocity relative to the atmosphere [km/s]

% Drag acceleration, 1e3 converting (kg/m^3)(m^2/kg)(km/s)^2 to km/s^2:
aD = -0.5e3 * rho * cfg.drag.Cd * cfg.drag.area / cfg.drag.mass * norm(vRel) * vRel;

end


function aS = accelSRP(t, r, cfg)
% accelSRP: perturbing acceleration from solar radiation pressure on a sphere, the spherical model of STK.
% t: time elapsed [s]
% r: position vector in ECI [km]
% cfg: parameter struct from config.m
% aS: solar radiation pressure acceleration in ECI [km/s^2]

rsun = sunVector(t, cfg); % rsun: Sun position in ECI [km]
d = rsun - r; % d: spacecraft to Sun vector [km]
dn = norm(d); % dn: spacecraft to Sun distance [km]
sh = rsun/norm(rsun); % sh: unit vector from the Earth centre to the Sun
along = dot(r, sh); % along: component of r along the Sun line [km]

% Earth shadow from shadow.m only behind the Earth and inside the penumbra cone, whose radius grows from cfg.Re by
% (cfg.Rsun + cfg.Re)/AU, under 0.005 km per km behind the Earth; everywhere else the whole Sun disk shows:
if along < 0 && norm(r - along*sh) < 1.01*cfg.Re - 0.01*along
    K = shadow(r, rsun, cfg); % K: visible fraction of the solar disk, 0 in umbra and 1 in full sunlight
else
    K = 1;
end

% Spherical model, a = -Cr A K S0/(m c R_AU^2) s_hat with R_AU the distance in AU and s_hat the unit vector to the Sun,
% 1e-3 converting (W/m^2)(m^2)/(kg m/s) = m/s^2 to km/s^2:
aS = -1e-3 * cfg.srp.Cr * cfg.srp.area * K * cfg.S0 / (cfg.drag.mass * cfg.c * (dn/cfg.AU)^2) * d/dn;

end
