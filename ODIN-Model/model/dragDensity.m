function rho = dragDensity(t, r, cfg)
% dragDensity: atmospheric mass density at points in the model frame, from the J77 table or from NRLMSISE-00 as cfg.drag.model sets.
% t: 1 x N times since sim start [s]
% r: 3 x N positions in the model frame [km]
% cfg: parameter struct from config.m, with cfg.drag.logRho from dragTable.m for J77
% rho: 1 x N mass density [kg/m^3]

rn = vecnorm(r); % rn: 1 x N distance from the centre of the Earth [km]

% Height above the WGS84 ellipsoid, first order in flattening:
h = rn - cfg.Re*(1 - cfg.f*(r(3,:)./rn).^2); % h: 1 x N height above the ellipsoid [km]

switch cfg.drag.model
    case 'J77'
        rho = exp(cfg.drag.logRho(h));

    case 'NRLMSISE-00'
        th = earthRotation(t, cfg); % th: 1 x N Greenwich angle [rad]

        % Model frame to Earth fixed, rotation about z by the Greenwich angle:
        xe = cos(th).*r(1,:) + sin(th).*r(2,:); % xe: 1 x N Earth fixed x [km]
        ye = cos(th).*r(2,:) - sin(th).*r(1,:); % ye: 1 x N Earth fixed y [km]

        e2 = cfg.f*(2 - cfg.f); % e2: squared eccentricity of the ellipsoid
        bE = cfg.Re*(1 - cfg.f); % bE: polar radius of the ellipsoid [km]
        p = hypot(xe, ye); % p: 1 x N distance from the polar axis [km]
        u = atan2(cfg.Re*r(3,:), bE*p); % u: 1 x N parametric latitude [rad]

        % Geodetic latitude, Bowring's formula:
        lat = atan2d(r(3,:) + e2/(1 - e2)*bE*sin(u).^3, p - e2*cfg.Re*cos(u).^3); % lat: 1 x N geodetic latitude [deg]

        lon = atan2d(ye, xe); % lon: 1 x N longitude, east positive [deg]

        dn = datenum(cfg.epoch) + t(:)/86400; % dn: N x 1 date number [days]
        d0 = floor(dn); % d0: N x 1 date number at 00:00 UTC of the day [days]
        Y = datevec(d0); % Y: N x 6 date vector at 00:00 UTC of the day
        doy = d0 - datenum(Y(:,1), 1, 0); % doy: N x 1 day of the year, 1 on 1 January
        sec = 86400*(dn - d0); % sec: N x 1 time of day, UTC [s]

        % Local solar time, UTC hours plus longitude over 15 deg/h:
        lst = sec/3600 + lon(:)/15; % lst: N x 1 local solar time [h]

        n = numel(h); % n: number of points
        [~, d] = atmosnrlmsise00(1e3*h(:), lat(:), lon(:), Y(:,1), doy, sec, lst, cfg.drag.F107, cfg.drag.F107, ...
                                 cfg.drag.Ap*ones(n, 7), ones(1, 23), 'Oxygen', 'None');
        % d: N x 9 densities from NRLMSISE-00, column 6 the total mass density with anomalous oxygen [kg/m^3]

        rho = d(:,6).';

    otherwise
        error('dragDensity:model', 'cfg.drag.model must be ''J77'' or ''NRLMSISE-00'', not ''%s''', cfg.drag.model);
end

end
