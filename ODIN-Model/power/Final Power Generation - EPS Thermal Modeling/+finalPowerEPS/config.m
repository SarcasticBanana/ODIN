function cfg = config()
% config: returns the parameter struct used by every other function.
% cfg: struct holding all constants, orbit elements and solver settings.

cfg=finalPowerEPS.modelSettings('orbit');
end
