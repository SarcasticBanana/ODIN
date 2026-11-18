function report=validateFinalPowerIntegration(repositoryRoot)
% VALIDATEFINALPOWERINTEGRATION Temporary native ODIN/package test environment.
% Tests run on local reference files only; cleanup restores path and folder.
oldPath=path; oldFolder=pwd;
cleanup=onCleanup(@()restore(oldPath,oldFolder)); %#ok<NASGU>
packageRoot=fileparts(fileparts(mfilename('fullpath')));
cd(repositoryRoot); addpath(repositoryRoot); odinPath;
addpath(packageRoot,fullfile(packageRoot,'tests'));
fprintf('Unchanged repository root: %s\n',repositoryRoot);
report.fixture=testFinalPowerEPS();
cfg=config();
fprintf('Workbook area %.6g m2; drag=%d, SRP=%d\n',cfg.power.panelArea,cfg.use.drag,cfg.use.srp);
out=propagate(cfg,1);
report.workbook=testFinalPowerEPS(cfg,out);
options=struct('makePlots',false,'showSummary',false,'refineOrbit',true, ...
    'arrayOverrides',struct('warnTemperatureExtrapolation',false,'warnLumpedModel',false,'warnProductTemperature',false));
originalOut=out; originalCfg=cfg; originalResolution=which('propagate');
r=runFinalPowerEPS(cfg,out,options);
assert(isequaln(out,originalOut) && isequaln(cfg,originalCfg));
assert(strcmp(which('propagate'),originalResolution));
assert(r.orbit.t(end)==out.t(end));
[isOriginal,index]=ismember(out.t,r.orbit.t);
assert(all(isOriginal));
assert(isequal(r.orbit.r(:,index),out.r) && isequal(r.orbit.v(:,index),out.v));
assert(all(r.orbit.illum>=0 & r.orbit.illum<=1));
assert(max(abs(r.thermal.nominal.energyResidual_W))<1e-7);
report.refinement=r.refinement;
report.effectiveArea_m2=r.mapping.effectiveArea_m2;
report.effectiveAgingRetained=r.mapping.effectiveAgingRetained;
report.nominalEnergyResidual_W=max(abs(r.thermal.nominal.energyResidual_W));
report.peakEnergyResidual_W=max(abs(r.thermal.peak.energyResidual_W));
report.complete=true;
fprintf('PASS: actual workbook scenario, original drag/SRP orbit, private refinement, input preservation and energy balance.\n');
end
function restore(oldPath,oldFolder)
path(oldPath); cd(oldFolder);
end
