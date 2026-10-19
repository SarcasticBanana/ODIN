function [sep, aSun, aEarth] = shadowGeometry(r, rsun, cfg)
% shadowGeometry: apparent disk radii and separation used by every shadow test.
% r: 3 x N satellite position in ECI [km]
% rsun: 3 x N Sun position in ECI [km]
% cfg: parameter struct from config.m
% sep: 1 x N angular separation of the Sun and Earth disk centres [rad]
% aSun: 1 x N apparent angular radius of the Sun [rad]
% aEarth: 1 x N apparent angular radius of the Earth [rad]

rs = rsun - r; % rs: 3 x N satellite to Sun vector [km]
rn = vecnorm(r); % rn: 1 x N orbital radius [km]
rsn = vecnorm(rs); % rsn: 1 x N satellite to Sun distance [km]

aSun = asin(min(cfg.Rsun./rsn, 1));
aEarth = asin(min(cfg.Re./rn, 1));

% Angle at the satellite between the Earth centre and the Sun centre:
sep = acos(min(max(dot(-r, rs, 1)./(rn.*rsn), -1), 1));

end
