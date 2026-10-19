function illum = shadow(r, rsun, cfg)
% shadow: fraction of the solar disk left visible by the Earth shadow.
% r: 3 x N satellite position in ECI [km]
% rsun: 3 x N Sun position in ECI [km]
% cfg: parameter struct from config.m
% illum: 1 x N visible fraction of the solar disk, 0 in umbra and 1 in full sunlight

switch lower(cfg.eclipse.model)
    case 'cylindrical'
        illum = cylindricalShadow(r, rsun, cfg);
    case 'conical'
        illum = conicalShadow(r, rsun, cfg);
    otherwise
        error('shadow:badModel','cfg.eclipse.model must be cylindrical or conical');
end

end


function illum = cylindricalShadow(r, rsun, cfg)
% cylindricalShadow: umbra only, a half infinite cylinder of radius cfg.Re behind the Earth.
% r: 3 x N satellite position in ECI [km]
% rsun: 3 x N Sun position in ECI [km]
% cfg: parameter struct from config.m
% illum: 1 x N visible fraction of the solar disk, 0 or 1

shat = rsun ./ vecnorm(rsun); % shat: 3 x N unit vector toward the Sun

along = dot(r, shat, 1); % along: 1 x N component of r along the Sun line [km]
across = vecnorm(r - along.*shat); % across: 1 x N offset of r from the Sun line [km]

illum = double(~(along < 0 & across < cfg.Re));

end


function illum = conicalShadow(r, rsun, cfg)
% conicalShadow: umbra and penumbra from the apparent disks of the Sun and the Earth.
% r: 3 x N satellite position in ECI [km]
% rsun: 3 x N Sun position in ECI [km]
% cfg: parameter struct from config.m
% illum: 1 x N visible fraction of the solar disk, 0 to 1

[sep, aSun, aEarth] = shadowGeometry(r, rsun, cfg);
% sep: 1 x N disk centre separation [rad]
% aSun: 1 x N apparent angular radius of the Sun [rad]
% aEarth: 1 x N apparent angular radius of the Earth [rad]

illum = ones(size(sep));

% Umbra, the Earth disk completely covers the Sun disk:
illum(sep <= aEarth - aSun) = 0;

% Antumbra, the Sun disk completely surrounds the Earth disk:
ant = sep <= aSun - aEarth; % ant: 1 x N mask of the annular eclipse samples
illum(ant) = 1 - (aEarth(ant)./aSun(ant)).^2;

% Penumbra, the overlap area of two circles of radii aSun and aEarth separated by sep:
pen = sep < aSun + aEarth & sep > abs(aEarth - aSun); % pen: 1 x N mask of the partial eclipse samples

x = (sep(pen).^2 + aSun(pen).^2 - aEarth(pen).^2) ./ (2*sep(pen));
% x: 1 x M distance from the Sun disk centre to the overlap chord [rad]
y = sqrt(max(aSun(pen).^2 - x.^2, 0)); % y: 1 x M half length of the overlap chord [rad]

cSun = min(max(x./aSun(pen), -1), 1); % cSun: 1 x M cosine argument of the Sun disk circular segment
cEarth = min(max((sep(pen)-x)./aEarth(pen), -1), 1); % cEarth: 1 x M cosine argument of the Earth disk circular segment

A = aSun(pen).^2.*acos(cSun) + aEarth(pen).^2.*acos(cEarth) - sep(pen).*y;
% A: 1 x M overlapped area of the two apparent disks [rad^2]

illum(pen) = 1 - A ./ (pi*aSun(pen).^2);

end
