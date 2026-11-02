function pm = powerModes(out, cfg, tilt)
% powerModes: battery state of charge and operating mode through the run, starting full in Nominal, the panel temperature carried step by step, loads and thresholds from the Operating modes section of the sheet.
% out: state history struct from propagate.m, cfg.use.eclipse on
% cfg: parameter struct from config.m
% tilt: 1 x N tilt of the panel normal from the Sun over the step from each output time, or one value for the whole run, cfg.slew.tiltSun when omitted [deg]
% pm: struct of 1 x N histories soc, mode (1 Nominal, 2 Raven, 3 Safe), pArr [W], pLoad [W] and T [K], and tMode, time in each mode [s]

N = numel(out.t); % N: number of output samples
if nargin < 3
    tilt = cfg.slew.tiltSun;
end
if isscalar(tilt)
    tilt = tilt*ones(1, N);
end

r = arrayPower(struct('t', out.t(1:2), 'illum', out.illum(1:2), 'rsun', out.rsun(:, 1:2)), cfg, cfg.power.panelArea, ...
               struct('incidenceAngles_deg', 0, 'incidenceSummaryAngles_deg', 0, 'warnTemperatureExtrapolation', false, 'warnLumpedModel', false));
% r: arrayPower over the first output step, run for its parameter struct
p = r.params; % p: arrayPower parameters, end-of-life factor from the sheet included

S = solarFlux(out, cfg); % S: 1 x N sunlight on a surface held normal to the Sun [W/m^2]

pm.soc = cfg.batt.socTop*ones(1, N); % pm.soc: 1 x N state of charge, share of rated capacity
pm.mode = ones(1, N); % pm.mode: 1 x N operating mode, 1 Nominal, 2 Raven, 3 Safe
pm.pArr = zeros(1, N); % pm.pArr: 1 x N array power at the bus, at the tilt of the step from each sample [W]
pm.T = (p.initialTemperature_C + 273.15)*ones(1, N); % pm.T: 1 x N panel temperature, arrayPower's starting value first [K]

for k = 1:N-1
    % Sunlight on the front face, the step's tilt held across it:
    f = S(k:k+1)*max(cosd(tilt(k)), 0); % f: sunlight on the front face at the start and end of the step [W/m^2]

    [pm.T(k+1), pS, pE] = panelStep(pm.T(k), f, out.t(k+1) - out.t(k), p, cfg.power.panelArea);
    % pS: array power at the bus at the start of the step [W]
    % pE: array power at the bus at the end of the step [W]
    pm.pArr(k) = pS;
    pm.pArr(k+1) = pE;

    % Array power at its mean over the step:
    [pm.soc(k+1), pm.mode(k+1)] = modeStep(pm.soc(k), pm.mode(k), (pS + pE)/2, out.t(k+1) - out.t(k), cfg);
end

pL = [cfg.modes.pNominal, cfg.modes.pRaven, cfg.modes.pSafe]; % pL: bus load of each mode [W]
pm.pLoad = pL(pm.mode); % pm.pLoad: 1 x N bus load [W]
pm.tMode = accumarray(pm.mode(1:N-1).', diff(out.t).', [3 1]).'; % pm.tMode: 1 x 3 time in Nominal, Raven and Safe [s]

end
