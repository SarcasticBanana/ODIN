function [logRho, Tinf] = dragTable(cfg)
% dragTable: J77 mass density against height for the exospheric temperature set by cfg.drag.F107.
% cfg: parameter struct from config.m
% logRho: pchip interpolant of the natural log of mass density [kg/m^3] against height above the ellipsoid [km]
% Tinf: exospheric temperature [K]

if ~exist('j77sri', 'file')
    addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'external', 'j77'));
end

F = cfg.drag.F107; % F: 10.7 cm solar flux, daily value and 81-day mean [sfu]

% J77 exospheric temperature at Kp = 0, the mean of the global diurnal extremes, daily flux equal to the 81-day mean:
Tinf = 5.48*F^0.8 + 101.8*F^0.4;

zMax = 2500; % zMax: top of the J77 profile [km]
[~, ~, ~, ~, ~, ~, ~, n, M] = j77sri(zMax, Tinf, 0:zMax);
% n: 1 x zMax+1 total number density at 0, 1, ..., zMax km [1/cm^3]
% M: 1 x zMax+1 mean molecular weight [g/mol]

NA = 6.02214076e23; % NA: Avogadro constant [1/mol]
h = 90:zMax; % h: J77 heights at the profile's 1 km spacing [km]

% Mass density from number density and mean molecular weight, 1e3 converting g/cm^3 to kg/m^3:
rho = n(h+1) .* M(h+1) * 1e3 / NA; % rho: J77 mass density at h [kg/m^3]

logRho = griddedInterpolant(h, log(rho), 'pchip', 'none');

end
