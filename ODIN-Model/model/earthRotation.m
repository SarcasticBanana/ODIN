function th = earthRotation(t, cfg)
% earthRotation: angle of the Greenwich meridian east of the x axis of the model frame, the mean equinox of cfg.epoch.
% t: 1 x N times since sim start [s]
% cfg: parameter struct from config.m
% th: 1 x N Greenwich angle [rad]

% Julian date from the serial date number:
jd0 = datenum(cfg.epoch) + 1721058.5; % jd0: Julian date at cfg.epoch, UTC taken as UT1 [days]

% Days and Julian centuries since the J2000.0 epoch:
d = jd0 - 2451545.0 + t(:)'/86400; % d: 1 x N time since J2000.0 [days]
T = d/36525; % T: 1 x N time since J2000.0 [Julian centuries]
T0 = (jd0 - 2451545.0)/36525; % T0: time from J2000.0 to cfg.epoch [Julian centuries]

% Greenwich mean sidereal time from the mean equinox of date, IAU 1982:
gmst = 280.46061837 + 360.98564736629*d + 0.000387933*T.^2 - T.^3/38710000; % gmst: 1 x N Greenwich mean sidereal time [deg]

% Mean equinox of date to that of cfg.epoch, general precession in right ascension of 1.2812 deg per century:
th = deg2rad(mod(gmst - 1.2812*(T - T0), 360));

end
