function dec = dragDecay(cfg, nDay)
% dragDecay: drag decay of a circular orbit, the semi-major axis stepped by the midpoint rule with the drag averaged around one orbit at each evaluation, so the day and night sides of the atmosphere both count.
% cfg: parameter struct from config.m, density source in cfg.drag.model, step, samples per orbit and floor in cfg.decay
% nDay: run length [days]
% dec: struct of 1 x N histories t [s], h [km], raan [deg], rho [kg/m^3] and dadt [km/s], and tFloor [s], the time h first falls below cfg.decay.hFloor, NaN if never

if strcmp(cfg.drag.model, 'J77')
    cfg.drag.logRho = dragTable(cfg);
    % cfg.drag.logRho: pchip interpolant of the natural log of J77 mass density [kg/m^3] against height above the ellipsoid [km]
end

dt = 86400*cfg.decay.dt; % dt: step length [s]
nStep = ceil(nDay*86400/dt); % nStep: number of steps

dec.t = (0:nStep)*dt; % dec.t: 1 x N times since sim start [s]
dec.h = NaN(1, nStep + 1); % dec.h: 1 x N semi-major axis less cfg.Re [km]
dec.raan = NaN(1, nStep + 1); % dec.raan: 1 x N right ascension of the ascending node [deg]
dec.rho = NaN(1, nStep + 1); % dec.rho: 1 x N mean density around the orbit [kg/m^3]
dec.dadt = NaN(1, nStep + 1); % dec.dadt: 1 x N rate of change of semi-major axis [km/s]
dec.tFloor = NaN; % dec.tFloor: time h first falls below cfg.decay.hFloor, linear between steps [s]

a = cfg.Re + cfg.h; % a: semi-major axis [km]
W = cfg.raan0; % W: right ascension of the ascending node [deg]
k = 1; % k: index of the current time
while true
    [dadt, dW, rhoBar] = orbitRates(dec.t(k), a, W, cfg);
    % dadt: orbit averaged rate of change of semi-major axis [km/s]
    % dW: rate of change of the node [deg/s]
    % rhoBar: mean density around the orbit [kg/m^3]

    dec.h(k) = a - cfg.Re;
    dec.raan(k) = mod(W, 360);
    dec.rho(k) = rhoBar;
    dec.dadt(k) = dadt;

    if dec.h(k) < cfg.decay.hFloor
        dec.tFloor = 0;
        if k > 1
            % Crossing of the floor, linear between the last two steps:
            dec.tFloor = dec.t(k-1) + (cfg.decay.hFloor - dec.h(k-1))/(dec.h(k) - dec.h(k-1))*dt;
        end
        break
    end
    if k > nStep
        break
    end

    % Midpoint rule:
    [dadtM, dWM] = orbitRates(dec.t(k) + dt/2, a + dt/2*dadt, W + dt/2*dW, cfg);
    % dadtM: rate of change of semi-major axis at the midpoint [km/s]
    % dWM: rate of change of the node at the midpoint [deg/s]

    a = a + dt*dadtM;
    W = W + dt*dWM;
    k = k + 1;
end

dec.t = dec.t(1:k);
dec.h = dec.h(1:k);
dec.raan = dec.raan(1:k);
dec.rho = dec.rho(1:k);
dec.dadt = dec.dadt(1:k);

end


function [dadt, dW, rhoBar] = orbitRates(t, a, W, cfg)
% orbitRates: drag rate of change of semi-major axis averaged around one circular orbit from time t, and the J2 rate of the node.
% t: time of the first sample since sim start [s]
% a: semi-major axis [km]
% W: right ascension of the ascending node at t [deg]
% cfg: parameter struct from config.m
% dadt: orbit averaged rate of change of semi-major axis [km/s]
% dW: rate of change of the node [deg/s]
% rhoBar: mean density around the orbit [kg/m^3]

n = sqrt(cfg.mu/a^3); % n: mean motion [rad/s]
inc = deg2rad(cfg.inc); % inc: inclination [rad]
k2 = 1.5*cfg.J2*(cfg.Re/a)^2; % k2: factor of the J2 secular rates

% Secular J2 rates of a circular orbit, node -k2 n cos i, argument of latitude n (1 + k2 (4 cos^2 i - 1)):
dW = -rad2deg(k2*n*cos(inc));
du = n*(1 + k2*(4*cos(inc)^2 - 1)); % du: rate of the argument of latitude [rad/s]

M = cfg.decay.nSample; % M: samples around the orbit
tau = (0:M-1)*2*pi/du/M; % tau: 1 x M sample times after t [s]

% Argument of latitude carried at the J2 rate from cfg.argp0 + cfg.nu0, node at the J2 rate within the orbit:
u = deg2rad(cfg.argp0 + cfg.nu0) + du*(t + tau); % u: 1 x M argument of latitude [rad]
Om = deg2rad(W + dW*tau); % Om: 1 x M right ascension of the ascending node [rad]

% Circular orbit, position and velocity from the node, inclination and argument of latitude:
r = a*[cos(Om).*cos(u) - sin(Om).*sin(u)*cos(inc)
       sin(Om).*cos(u) + cos(Om).*sin(u)*cos(inc)
       sin(u)*sin(inc)]; % r: 3 x M positions in the model frame [km]
v = sqrt(cfg.mu/a)*[-cos(Om).*sin(u) - sin(Om).*cos(u)*cos(inc)
                    -sin(Om).*sin(u) + cos(Om).*cos(u)*cos(inc)
                    cos(u)*sin(inc)]; % v: 3 x M velocities in the model frame [km/s]

rho = dragDensity(t + tau, r, cfg); % rho: 1 x M mass density [kg/m^3]

% Velocity relative to an atmosphere rotating rigidly with the Earth, v - wE x r:
vRel = v - cfg.wE*[-r(2,:); r(1,:); zeros(1, M)]; % vRel: 3 x M velocity relative to the atmosphere [km/s]

% Drag acceleration as in eom.m, 1e3 converting (kg/m^3)(m^2/kg)(km/s)^2 to km/s^2:
aD = -0.5e3*cfg.drag.Cd*cfg.drag.area/cfg.drag.mass*rho.*vecnorm(vRel).*vRel; % aD: 3 x M drag acceleration [km/s^2]

% Semi-major axis rate from the specific energy rate, da/dt = (2 a^2/mu) <v . aD>:
dadt = 2*a^2/cfg.mu*mean(sum(v.*aD, 1));

rhoBar = mean(rho);

end
