function R = gcrfRotation(cfg)
% gcrfRotation: rotation from the mean equator and equinox of cfg.epoch to GCRF, IAU 2006 precession with frame bias.
% cfg: parameter struct from config.m
% R: 3x3 rotation matrix, r_GCRF = R * r

dTT = 69.184; % dTT: TT minus UTC, 37 leap seconds plus 32.184 s [s]
as = pi/648000; % as: one arcsecond [rad]

% Julian centuries of TT since J2000.0:
T = (juliandate(datetime(cfg.epoch)) + dTT/86400 - 2451545.0) / 36525; % T: time argument of the precession series [Julian centuries]

% IAU 2006 Fukushima-Williams precession angles, frame bias included:
gamb = as * polyval([ 0.0000000260 -0.000002788 -0.00031238  0.4932044   10.556378    -0.052928], T); % gamb: GCRS right ascension of the ecliptic of date node [rad]
phib = as * polyval([-0.0000000176 -0.000000440  0.00053289  0.0511268  -46.811016 84381.412819], T); % phib: inclination of the ecliptic of date to the GCRS equator [rad]
psib = as * polyval([-0.0000000148 -0.000026452 -0.00018522  1.5584175 5038.481484    -0.041775], T); % psib: arc of the ecliptic of date from that node to the mean equinox [rad]
epsa = as * polyval([-0.0000000434 -0.000000576  0.00200340 -0.0001831  -46.836769 84381.406   ], T); % epsa: mean obliquity of date [rad]

% GCRS to mean of date bias precession matrix, transposed:
R = (R1(-epsa) * R3(-psib) * R1(phib) * R3(gamb)).';

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
