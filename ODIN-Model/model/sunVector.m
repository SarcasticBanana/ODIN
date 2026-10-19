function [rsun, sunhat] = sunVector(t, cfg)
% sunVector: low precision Sun position in ECI over the propagation window.
% t: 1 x N times since sim start [s]
% cfg: parameter struct from config.m
% rsun: 3 x N Sun position vector in ECI [km]
% sunhat: 3 x N unit vector from the Earth centre toward the Sun

jd = julianDate(cfg.epoch) + t(:)'/86400; % jd: 1 x N Julian date at each output time [days]

% Julian centuries elapsed since the J2000.0 epoch:
T = (jd - 2451545.0) / 36525; % T: 1 x N time argument of the series [Julian centuries]

lamM = deg2rad(280.460 + 36000.771*T); % lamM: 1 x N mean ecliptic longitude of the Sun [rad]
M = deg2rad(357.5291092 + 35999.05034*T); % M: 1 x N mean anomaly of the Sun [rad]

% Equation of the centre added to the mean longitude:
lam = lamM + deg2rad(1.914666471)*sin(M) + deg2rad(0.019994643)*sin(2*M);
% lam: 1 x N apparent ecliptic longitude of the Sun [rad]

obl = deg2rad(23.439291 - 0.0130042*T); % obl: 1 x N mean obliquity of the ecliptic [rad]

% Earth to Sun distance from the low order eccentricity expansion:
rmag = (1.000140612 - 0.016708617*cos(M) - 0.000139589*cos(2*M)) * cfg.AU;
% rmag: 1 x N Earth to Sun distance [km]

% Ecliptic longitude rotated into equatorial components about the x axis:
sunhat = [cos(lam)
          cos(obl).*sin(lam)
          sin(obl).*sin(lam)];

rsun = sunhat .* rmag;

end


function jd = julianDate(dv)
% julianDate: Julian date of a UTC calendar instant.
% dv: [year month day hour minute second] in UTC
% jd: Julian date [days]

y = dv(1); % y: calendar year
m = dv(2); % m: calendar month
d = dv(3); % d: calendar day of the month

if m <= 2
    y = y - 1;
    m = m + 12;
end

c = floor(y/100); % c: century number of the shifted calendar year
b = 2 - c + floor(c/4); % b: Gregorian calendar correction [days]

dayFrac = (dv(4) + dv(5)/60 + dv(6)/3600) / 24; % dayFrac: elapsed fraction of the calendar day

% Gregorian calendar date to Julian date:
jd = floor(365.25*(y + 4716)) + floor(30.6001*(m + 1)) + d + b - 1524.5 + dayFrac;

end
