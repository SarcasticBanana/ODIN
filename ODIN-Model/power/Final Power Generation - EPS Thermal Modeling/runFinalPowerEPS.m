function result=runFinalPowerEPS(cfg,out,options)
% RUNFINALPOWEREPS Final Power Generation / EPS Thermal Modeling.
% In ODIN: odinPath; cfg=config(); out=propagate(cfg,10);
% Add ONLY this folder (the parent of +finalPowerEPS), then:
% result=runFinalPowerEPS(cfg,out);
% No shared functions are replaced, and no path/current-folder changes occur.
if nargin<2
    error('finalPowerEPS:inputs','Supply cfg and out from your ODIN run. See the package README for setup.');
end
if nargin<3, options=struct(); end
result=finalPowerEPS.run(cfg,out,options);
end
