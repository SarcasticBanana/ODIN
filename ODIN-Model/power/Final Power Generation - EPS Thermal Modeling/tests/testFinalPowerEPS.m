function report=testFinalPowerEPS(cfg,out)
% TESTFINALPOWEREPS Regression and isolation checks against the unchanged ODIN code.
% First run odinPath in the ODIN repository, then add the parent of +finalPowerEPS
% and this tests folder to the MATLAB path. This test never changes path or pwd.
% No arguments uses an explicit test-only orbit fixture (drag/SRP off). To
% check the actual workbook scenario, pass cfg=config() and out=propagate(cfg,1).
% These fixture values do not change any production model or workbook defaults.

assert(~isempty(which('finalPowerEPS.run')),'Add the package parent folder before running this test.');
teamRoot=fileparts(which('odinPath'));
assert(~isempty(teamRoot),'Run odinPath from the intended unchanged ODIN checkout first.');
shared={'config','propagate','eom','shadow','shadowGeometry','eclipseTimes', ...
    'sunVector','solarFlux','oe2rv','arrayPower','panelStep','modeStep','powerModes','slewFrequency'};
beforeResolution=cellfun(@which,shared,'UniformOutput',false);
for k=1:numel(shared)
    assert(~isempty(beforeResolution{k}),['Missing team function: ' shared{k}]);
    assert(startsWith(lower(beforeResolution{k}),lower([teamRoot filesep])), ...
        ['Team function is shadowed before the test: ' shared{k} ' -> ' beforeResolution{k}]);
end
originalPath=path; originalFolder=pwd;
if nargin==0
    cfg=fixture();
end
if nargin<2
    out=propagate(cfg,1);
end
cfgBefore=cfg; outBefore=out;
figureCount=numel(findall(groot,'Type','figure'));
legacyBefore=legacyFingerprint(out,cfg);

arrayOverrides=struct('runIncidenceSweep',false,'makePlots',false, ...
    'makeMainArrayPlots',false,'showAreaSummary',false,'showIncidenceSummary',false, ...
    'warnTemperatureExtrapolation',false,'warnLumpedModel',false,'warnProductTemperature',false);
options=struct('makePlots',false,'showSummary',false,'refineOrbit',false,'includePowerBalance',false, ...
    'arrayOverrides',arrayOverrides,'componentOverrides',struct('showSummary',false));
result=finalPowerEPS.run(cfg,out,options);
count=0;
check(isequaln(cfg,cfgBefore),'Input cfg was changed.');
check(isequaln(out,outBefore),'Input orbit was changed.');
check(strcmp(path,originalPath),'MATLAB path was changed.');
check(strcmp(pwd,originalFolder),'MATLAB current folder was changed.');
afterResolution=cellfun(@which,shared,'UniformOutput',false);
check(isequal(beforeResolution,afterResolution),'A shared team function changed resolution.');
check(numel(findall(groot,'Type','figure'))==figureCount,'makePlots=false created a figure.');
check(isequal(result.array.time_s,out.t),'refineOrbit=false changed array sampling.');
check(isequal(result.orbit.t,out.t),'refineOrbit=false changed the returned orbit time grid.');
check(result.array.params.thermalNodeCount==2,'Updated analysis lost its two thermal nodes.');
check(result.array.params.backEmissivity>0,'Rear-face radiation was disabled.');
check(max(result.array.thermal.backRadiation_W_m2)>0,'Rear-face emission was not computed.');
check(isfield(result,'mapping') && ~isempty(result.mapping),'Configuration provenance is missing.');
if isfield(cfg,'power') && isfield(cfg.power,'panelArea')
    check(abs(result.array.params.selectedArrayArea_m2-cfg.power.panelArea)<1e-12, ...
        'Workbook area was not mapped to the updated analysis.');
end
eolFields={'eolRadiation','eolUV','eolDebris','eolThermalCycling','eolContamination'};
if isfield(cfg,'power') && all(isfield(cfg.power,eolFields))
    eol=1;
    for k=1:numel(eolFields), eol=eol*cfg.power.(eolFields{k}); end
    check(abs(result.array.params.agingRetained-eol)<1e-12, ...
        'The workbook EOL factors were not mapped consistently.');
end
for mode={'nominal','peak'}
    label=mode{1}; v=result.thermal.(label); iv=result.iv.(label);
    check(max(abs(v.energyResidual_W))<1e-7,[label ' total energy does not close.']);
    check(max(abs(v.solar.balanceResidual_W))<1e-7,[label ' thermal energy does not close.']);
    check(all(v.solar.electricalPower_W<=v.solar.potentialElectricalPower_W+1e-7), ...
        [label ' extraction exceeds available panel power.']);
    check(max(abs(iv.Pmp_W-v.solar.potentialElectricalPower_W))<1e-7, ...
        [label ' I-V proxy does not match available panel-terminal power.']);
    check(all(isfinite(iv.Vmp_V)) && all(iv.Vmp_V>=0),[label ' invalid I-V voltage.']);
end

legacyAfter=legacyFingerprint(out,cfg);
legacyNames=fieldnames(legacyBefore);
for k=1:numel(legacyNames)
    name=legacyNames{k}; a=legacyBefore.(name); b=legacyAfter.(name);
    check(isequal(size(a),size(b)) && max(abs(a(:)-b(:)))<1e-10, ...
        ['Unchanged ODIN calculation changed: ' name]);
end
check(isequal(beforeResolution,cellfun(@which,shared,'UniformOutput',false)), ...
    'Original function resolution changed during regression checks.');
check(strcmp(path,originalPath) && strcmp(pwd,originalFolder), ...
    'Regression checks changed path or current folder.');
report=struct('checksPassed',count,'teamRoot',teamRoot, ...
    'inputDuration_s',out.t(end)-out.t(1),'sharedFunctions',{shared}, ...
    'sharedResolution',{afterResolution},'originalCalculationsUnchanged',true);
fprintf('PASS: %d package integration, conservation and unchanged-team-code checks.\n',count);
    function check(ok,message)
        assert(all(ok(:)),message); count=count+1;
    end
end

function r=legacyFingerprint(out,cfg)
% Exercise actual call shapes used by the repository's powerModes/MISSION code.
subset=struct('t',out.t(1:2),'illum',out.illum(1:2),'rsun',out.rsun(:,1:2));
quiet=struct('incidenceAngles_deg',0,'incidenceSummaryAngles_deg',0, ...
    'warnTemperatureExtrapolation',false,'warnLumpedModel',false, ...
    'makePlots',false,'makeMainArrayPlots',false);
a=arrayPower(subset,cfg,cfg.power.panelArea,quiet);
r.arrayTemperature=a.cellTemperature_C;
r.arrayBusPower=a.busPlanning_W;
[r.panelTemperature,r.panelStartPower,r.panelEndPower]= ...
    panelStep(300,[1000 800],60,a.params,cfg.power.panelArea);
[r.modeSOC,r.modeValue]=modeStep(.8,1,20,60,cfg);
m=powerModes(out,cfg,0);
r.powerModeSOC=m.soc; r.powerMode=m.mode;
r.powerModeArray=m.pArr; r.powerModeTemperature=m.T;
end

function cfg=fixture()
% Explicit test case only. Native ODIN model functions perform propagation.
cfg=finalPowerEPS.modelSettings('orbit');
cfg.h=600; cfg.use.drag=false; cfg.use.srp=false;
cfg.ode.dtOut=60;
cfg.power=struct('panelArea',.175,'eolRadiation',.97,'eolUV',.99, ...
    'eolDebris',.98,'eolThermalCycling',1,'eolContamination',1,'years',5);
cfg.slew=struct('tiltSun',0,'tilt',63.3);
cfg.batt=struct('energy',45,'socTop',1,'etaChg',.9,'etaDis',.95);
cfg.modes=struct('pNominal',23.65,'pRaven',15,'pSafe',5, ...
    'socRaven',.6,'socSafe',.5,'socNominal',.95);
end
