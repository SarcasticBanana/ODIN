function result=run(cfg,out,options)
% RUN Final Power Generation / EPS Thermal Modeling, isolated ODIN extension.
% result=finalPowerEPS.run(cfg,out,options) consumes the team's configuration
% and orbit without modifying them, the MATLAB path, or the working folder.
% Existing RUN/MISSION/powerModes/panelStep continue using their original code.
if nargin<2
    error('finalPowerEPS:inputs','Supply cfg and out from your ODIN run. Example: cfg=config(); out=propagate(cfg,10); result=runFinalPowerEPS(cfg,out);');
end
if nargin<3, options=struct(); end
validateattributes(cfg,{'struct'},{'scalar'}); validateattributes(out,{'struct'},{'scalar'});
validateattributes(options,{'struct'},{'scalar'});
defaults=struct('makePlots',true,'refineOrbit',true,'includePowerBalance',false, ...
    'showSummary',true,'arrayOverrides',struct(),'componentOverrides',struct(), ...
    'ivOverrides',struct(),'powerBalanceOverrides',struct(),'refinementOptions',struct());
for field=fieldnames(options).'
    name=field{1};
    if ~isfield(defaults,name), error('finalPowerEPS:option','Unknown option: %s',name); end
    defaults.(name)=options.(name);
end
options=defaults;
for field={'makePlots','refineOrbit','includePowerBalance','showSummary'}
    value=options.(field{1});
    validateattributes(double(value),{'numeric'},{'real','finite','scalar','integer','>=',0,'<=',1});
end
for field={'arrayOverrides','componentOverrides','ivOverrides','powerBalanceOverrides','refinementOptions'}
    validateattributes(options.(field{1}),{'struct'},{'scalar'});
end
if isfield(options.ivOverrides,'mode') && ~strcmp(options.ivOverrides.mode,'power_consistent_proxy')
    error('finalPowerEPS:ivMode',['The coupled run uses power_consistent_proxy. ' ...
        'Run finalPowerEPS.solarArrayIV separately for the counterfactual xtj_reference study.']);
end
for field={'t','r','v','T'}
    if ~isfield(out,field{1}), error('finalPowerEPS:orbit','Orbit is missing out.%s.',field{1}); end
end
if options.refineOrbit
    [analysisOrbit,refinement]=finalPowerEPS.refineOrbit(out,cfg,options.refinementOptions);
else
    analysisOrbit=out;
    refinement=struct('enabled',false,'description','Original orbit grid retained explicitly.');
end
if ~isfield(analysisOrbit,'rsun') || ~isfield(analysisOrbit,'sunhat')
    [analysisOrbit.rsun,analysisOrbit.sunhat]=finalPowerEPS.sunVector(analysisOrbit.t,cfg);
end
if ~isfield(analysisOrbit,'illum')
    if isfield(cfg,'use') && isfield(cfg.use,'eclipse') && ~cfg.use.eclipse
        analysisOrbit.illum=ones(size(analysisOrbit.t));
    else
        analysisOrbit.illum=finalPowerEPS.shadow(analysisOrbit.r,analysisOrbit.rsun,cfg);
    end
end

% Map only the repository fields explicitly supported by this analysis.
% Different battery and mode models remain separate, with no silent conversion.
settings=finalPowerEPS.modelSettings('array');
mapping=struct('areaSource','package default','agingSource','package default', ...
    'originalPackageArea_m2',settings.selectedArrayArea_m2, ...
    'originalPackageAgingRetained',settings.agingRetained, ...
    'unmappedModels',{{'cfg.batt: original ODIN battery/aging model remains separate', ...
    'cfg.modes: original Nominal/Raven/Safe controller remains separate', ...
    'ground contacts: package component model retains its labeled periodic schedule'}});
changes=struct('runIncidenceSweep',false);
if isfield(cfg,'power') && isfield(cfg.power,'panelArea')
    validateattributes(cfg.power.panelArea,{'numeric'},{'real','finite','scalar','positive'});
    changes.selectedArrayArea_m2=cfg.power.panelArea;
    changes.incidenceSweepArea_m2=cfg.power.panelArea;
    mapping.areaSource='cfg.power.panelArea (repository workbook/configuration)';
end
eol={'eolRadiation','eolUV','eolDebris','eolThermalCycling','eolContamination'};
if isfield(cfg,'power') && any(isfield(cfg.power,eol))
    if ~all(isfield(cfg.power,eol))
        error('finalPowerEPS:eol','Supply all five repository EOL factors, or none; partial factors are ambiguous.');
    end
    factor=1;
    for k=1:numel(eol)
        value=cfg.power.(eol{k}); validateattributes(value,{'numeric'},{'real','finite','scalar','>=',0,'<=',1});
        factor=factor*value;
    end
    changes.agingRetained=factor;
    mapping.agingSource='Product of the five cfg.power EOL factors';
end
for field=fieldnames(options.arrayOverrides).'
    changes.(field{1})=options.arrayOverrides.(field{1});
end
if isfield(options.arrayOverrides,'selectedArrayArea_m2')
    mapping.areaSource='options.arrayOverrides.selectedArrayArea_m2';
    if ~isfield(options.arrayOverrides,'incidenceSweepArea_m2')
        changes.incidenceSweepArea_m2=changes.selectedArrayArea_m2;
    end
end
if isfield(options.arrayOverrides,'agingRetained'), mapping.agingSource='options.arrayOverrides.agingRetained'; end
% Explicit presentation switches apply to this invocation only.
changes.makePlots=logical(options.makePlots); changes.makeMainArrayPlots=false;
changes.showAreaSummary=false; changes.showIncidenceSummary=false;
array=finalPowerEPS.arrayPower(analysisOrbit,cfg,[],changes);
componentChanges=options.componentOverrides;
componentChanges.showSummary=logical(options.showSummary);
thermal=finalPowerEPS.componentThermalOrbit(analysisOrbit,cfg,array,componentChanges);
ivChanges=struct('arrayArea_m2',thermal.nominal.solar.params.area_m2, ...
    'cellCoverage',thermal.nominal.solar.params.cellCoverage);
for field=fieldnames(options.ivOverrides).'
    ivChanges.(field{1})=options.ivOverrides.(field{1});
end
ivChanges.makePlots=logical(options.makePlots); ivChanges.showSummary=logical(options.showSummary);
for mode={'nominal','peak'}
    label=mode{1}; solar=thermal.(label).solar;
    ivInput=struct('time_s',thermal.time_s,'irradiance_W_m2',solar.irradiance_W_m2, ...
        'cellTemperature_C',solar.cellTemperature_C,'availablePower_W',solar.potentialElectricalPower_W);
    iv.(label)=finalPowerEPS.solarArrayIV(ivInput,ivChanges);
    ivChanges.makePlots=false; ivChanges.showSummary=false;
end
result=struct('name','Final Power Generation / EPS Thermal Modeling', ...
    'configuration',cfg,'orbit',analysisOrbit,'refinement',refinement,'array',array, ...
    'thermal',thermal,'iv',iv,'mapping',mapping,'options',options,'figures',gobjects(0));
result.mapping.effectiveArea_m2=array.params.selectedArrayArea_m2;
result.mapping.effectiveAgingRetained=array.params.agingRetained;
result.mapping.assumedWiring=[iv.nominal.params.seriesCells iv.nominal.params.parallelStrings];
if options.includePowerBalance
    balanceOptions=options.powerBalanceOverrides;
    balanceArray=array.params;
    if isfield(balanceOptions,'arrayOverrides')
        validateattributes(balanceOptions.arrayOverrides,{'struct'},{'scalar'});
        for field=fieldnames(balanceOptions.arrayOverrides).'
            balanceArray.(field{1})=balanceOptions.arrayOverrides.(field{1});
        end
    end
    balanceOptions.arrayOverrides=balanceArray;
    balanceOptions.makePlot=logical(options.makePlots); balanceOptions.showSummary=logical(options.showSummary);
    result.powerBalance=finalPowerEPS.powerBalance(cfg,balanceOptions);
end
if options.makePlots
    result.figures=finalPowerEPS.plotModelDiagnostics(analysisOrbit,cfg,array,thermal,iv);
end
if options.showSummary
    fprintf('\n%s\nArea %.6g m2 from %s\nEOL factor %.6g from %s\n', ...
        result.name,result.mapping.effectiveArea_m2,result.mapping.areaSource, ...
        result.mapping.effectiveAgingRetained,result.mapping.agingSource);
    fprintf('Existing ODIN operations and battery-aging models remain unchanged.\n');
end
end
