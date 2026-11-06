function m = energyMargin(t, pArr, pLoad, tSh, cfg)
% energyMargin: energy margin of constant bus loads over one orbit, M_E = (E_gen - E_load)/E_load as on the metric definitions slide, the battery's net energy over the orbit and the spare power, the extra constant load that brings that net energy to zero.
% t: 1 x N times over one orbit [s]
% pArr: 1 x N array power at the bus [W]
% pLoad: 1 x M constant bus loads [W]
% tSh: shadow length inside the orbit [s]
% cfg: parameter struct, charge and discharge efficiencies from cfg.batt
% m: struct of T [s], eGen [Wh], and 1 x M eLoad [Wh], M, net [Wh] and spare [W]

m.T = t(end) - t(1); % m.T: orbit length [s]

% Energy generated over the orbit, trapezoidal:
m.eGen = trapz(t, pArr)/3600; % m.eGen: array energy at the bus over the orbit [Wh]

% Load energy with the shadow share through the battery, E_load = P (T - t_sh) + P t_sh/(eta_c eta_d):
m.eLoad = pLoad*((m.T - tSh) + tSh/(cfg.batt.etaChg*cfg.batt.etaDis))/3600; % m.eLoad: 1 x M load energy over the orbit [Wh]

m.M = (m.eGen - m.eLoad)./m.eLoad; % m.M: 1 x M energy margin

pm = (pArr(1:end-1) + pArr(2:end))/2; % pm: array power at the mean of each step [W]
dt = diff(t); % dt: step lengths [s]

% Net battery energy at a constant load P, surplus times the charge efficiency, deficit over the discharge efficiency:
net = @(P) sum(max(pm - P, 0)*cfg.batt.etaChg.*dt + min(pm - P, 0)/cfg.batt.etaDis.*dt)/3600; % net: battery energy gained over the orbit at load P [Wh]

m.net = arrayfun(net, pLoad); % m.net: 1 x M battery energy gained over the orbit [Wh]
pMax = fzero(net, [0 max(pArr)]); % pMax: constant load with no net battery energy over the orbit [W]
m.spare = pMax - pLoad; % m.spare: 1 x M extra constant load the orbit can carry [W]

end
