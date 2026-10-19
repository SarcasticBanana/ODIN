function [r, v] = oe2rv(a, ecc, incDeg, raanDeg, argpDeg, nuDeg, mu)
% oe2rv: converts classical orbital elements to ECI position and velocity.
% a: semi-major axis [km]
% ecc: eccentricity
% incDeg: inclination [deg]
% raanDeg: right ascension of the ascending node [deg]
% argpDeg: argument of perigee [deg]
% nuDeg: true anomaly [deg]
% mu: gravitational parameter [km^3/s^2]
% r: position vector in ECI [km]
% v: velocity vector in ECI [km/s]

i = deg2rad(incDeg); % i: inclination [rad]
Om = deg2rad(raanDeg); % Om: right ascension of the ascending node [rad]
w = deg2rad(argpDeg); % w: argument of perigee [rad]
nu = deg2rad(nuDeg); % nu: true anomaly [rad]

p = a * (1 - ecc^2); % p: semilatus rectum [km]

if p <= 0
    error('oe2rv:badOrbit','semi-latus rectum must be positive');
end

% Orbit equation in the perifocal frame:
rPQW = [p*cos(nu) / (1 + ecc*cos(nu))
        p*sin(nu) / (1 + ecc*cos(nu))
        0]; % rPQW: position in the perifocal frame [km]

% Perifocal velocity from the vis-viva and angular momentum relations:
vPQW = sqrt(mu/p) * [-sin(nu)
                     ecc + cos(nu)
                     0]; % vPQW: velocity in the perifocal frame [km/s]

% Perifocal to ECI rotation sequence:
R = R3(-Om) * R1(-i) * R3(-w); % R: perifocal to ECI direction cosine matrix
r = R * rPQW;
v = R * vPQW;

end


function M = R1(th)
% R1: frame rotation about the x axis.
% th: rotation angle [rad]
% M: 3x3 rotation matrix

M = [1 0 0
     0 cos(th) sin(th)
     0 -sin(th) cos(th)];

end


function M = R3(th)
% R3: frame rotation about the z axis.
% th: rotation angle [rad]
% M: 3x3 rotation matrix

M = [cos(th) sin(th) 0
     -sin(th) cos(th) 0
     0 0 1];

end
