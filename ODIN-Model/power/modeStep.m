function [soc, mode] = modeStep(soc, mode, pArr, dt, cfg)
% modeStep: battery state of charge and operating mode one step later, array power held over the step, loads and thresholds from the Operating modes section.
% soc: state of charge at the start of the step, share of full charge (the energy held at cfg.batt.socTop)
% mode: operating mode at the start of the step, 1 Nominal, 2 Raven, 3 Safe
% pArr: array power at the bus over the step [W]
% dt: step length [s]
% cfg: parameter struct from config.m
% soc: state of charge at the end of the step, share of full charge
% mode: operating mode at the end of the step

pLoad = [cfg.modes.pNominal, cfg.modes.pRaven, cfg.modes.pSafe]; % pLoad: bus load of each mode [W]
pBat = pArr - pLoad(mode); % pBat: array surplus, negative when the battery carries the load [W]

% Battery energy change, surplus times the charge efficiency, deficit over the discharge efficiency:
if pBat >= 0
    dE = pBat*cfg.batt.etaChg*dt/3600; % dE: battery energy change over the step [Wh]
else
    dE = pBat/cfg.batt.etaDis*dt/3600;
end

% Full charge, the rated energy up to the top of charge:
eFull = cfg.batt.energy*cfg.batt.socTop; % eFull: battery energy at full charge [Wh]

% Charge held between empty and full:
soc = min(max(soc + dE/eFull, 0), 1);

% Ladder down Nominal, Raven, Safe below socRaven and socSafe, up Safe to Raven above socRaven and Raven to Nominal above socNominal:
if mode == 1 && soc < cfg.modes.socRaven
    mode = 2;
end
if mode == 2 && soc < cfg.modes.socSafe
    mode = 3;
elseif mode == 2 && soc > cfg.modes.socNominal
    mode = 1;
elseif mode == 3 && soc > cfg.modes.socRaven
    mode = 2;
end

end
