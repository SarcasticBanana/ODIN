function [T, pS, pE] = panelStep(T, S, dt, p, area)
% panelStep: panel temperature and array power at the bus over one step, the plate and cell equations of arrayPower with the sunlight on the front face linear across the step.
% T: panel temperature at the start of the step [K]
% S: 1 x 2 sunlight on the front face at the start and end of the step [W/m^2]
% dt: step length [s]
% p: parameter struct of arrayPower, its result.params
% area: panel front-face area [m^2]
% T: panel temperature at the end of the step [K]
% pS: array power at the bus at the start of the step [W]
% pE: array power at the bus at the end of the step [W]

n = max(1, ceil(dt/p.thermalMaxStep_s)); % n: substeps, none longer than arrayPower's largest thermal step
h = dt/n; % h: substep length [s]
C = p.panelDensity_kg_m3*p.panelThickness_m*p.panelSpecificHeat_J_kgK; % C: areal heat capacity of the plate [J/(m^2 K)]

pS = bus(S(1), T, p, area);

for m = 1:n
    f0 = S(1) + (S(2) - S(1))*(m - 1)/n; % f0: sunlight at the start of the substep [W/m^2]
    fh = S(1) + (S(2) - S(1))*(m - 0.5)/n; % fh: sunlight at the middle of the substep [W/m^2]
    f1 = S(1) + (S(2) - S(1))*m/n; % f1: sunlight at the end of the substep [W/m^2]

    % Classical fourth-order Runge-Kutta on C dT/dt = absorbed - extracted - radiated:
    k1 = heat(T, f0, p)/C; % k1: temperature rate at the start [K/s]
    k2 = heat(T + h/2*k1, fh, p)/C; % k2: first midpoint rate [K/s]
    k3 = heat(T + h/2*k2, fh, p)/C; % k3: second midpoint rate [K/s]
    k4 = heat(T + h*k3, f1, p)/C; % k4: end rate [K/s]
    T = T + h/6*(k1 + 2*k2 + 2*k3 + k4);
end

pE = bus(S(2), T, p, area);

end


function q = heat(T, f, p)
% heat: net heat into the plate per unit face area, plateHeatRate of arrayPower.
% T: panel temperature [K]
% f: sunlight on the front face [W/m^2]
% p: parameter struct of arrayPower
% q: absorbed sunlight less the extracted electrical power and the front-face radiation [W/m^2]

sigma = 5.670374419e-8; % sigma: Stefan-Boltzmann constant [W/(m^2 K^4)]

% Absorbed sunlight, power taken out through the MPPT, radiation to the cold sink:
q = p.solarAbsorptivity*f*p.selfShadowRetained - density(f, T, p)*p.mpptEfficiency ...
    - p.frontEmissivity*sigma*(T^4 - p.spaceTemperature_K^4);

end


function a = density(f, T, p)
% density: electrical power per unit face area at the cells, arrayPowerDensity of arrayPower.
% f: sunlight on the front face [W/m^2]
% T: panel temperature [K]
% p: parameter struct of arrayPower
% a: available electrical power [W/m^2]

eff = p.cellEfficiency; % eff: cell efficiency at temperature T
if f > 0
    % Linear temperature correction about the reference temperature:
    eff = p.cellEfficiency*(1 + p.powerTempCoeff_perC*(T - 273.15 - p.referenceTemperature_C));
end

a = f*eff*p.cellCoverage*p.selfShadowRetained*p.mismatchRetained*p.opticalRetained;

end


function b = bus(f, T, p, area)
% bus: array power at the bus, end-of-life planning value as arrayPower's busPlanning_W.
% f: sunlight on the front face [W/m^2]
% T: panel temperature [K]
% p: parameter struct of arrayPower
% area: panel front-face area [m^2]
% b: array power at the bus [W]

b = area*density(f, T, p)*p.wiringDiodeRetained*p.mpptEfficiency*p.converterEfficiency*p.agingRetained;

end
