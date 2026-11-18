function result = SolarArrayThermalModel(inputs, overrides)
% SOLARARRAYTHERMALMODEL Total panel heat rejection, heat input and transient storage [W].
% Primary output: Qout_W = radiation to space/Earth + radiation/conduction to structure.
% QheatInput_W = absorbed solar/Earth - terminal electrical export + heater.
% Qstorage_W = total C*d(mean front/rear T)/dt; QheatInput_W = Qout_W + Qstorage_W.
% Resistive/conversion heating is already included; no separate percentage loss is added.

%% Usage
% s = finalPowerEPS.SolarArrayThermalModel(arrayResult); % postprocess arrayPower temperature history
% in = struct('timeEdges_s',0:10:6000,'irradiance_W_m2',S_interval);
% s = finalPowerEPS.SolarArrayThermalModel(in); % solve temperature with available electrical extraction
% in.electricalPower_W = actualPanelExport_W; % optional actual export BEFORE remote EPS losses
% s = finalPowerEPS.SolarArrayThermalModel(in); % solve with actual electrical extraction
% With timeEdges_s, outputs are interval averages; temperatureEdges_C gives boundary temperatures.
% Supplying cellTemperature_C instead postprocesses it; no new temperature is solved.
% time_s sample histories require supplied temperature. No arguments evaluates a reference point.

% Shared defaults and assumptions: edit finalPowerEPS.modelSettings('solar').
p=finalPowerEPS.modelSettings('solar');

%% Inherited settings and inputs
if nargin<1 || isempty(inputs), inputs=struct(); end
if nargin<2, overrides=struct(); end
checkStruct(inputs,'inputs'); checkStruct(overrides,'overrides');
source='supplied inputs or reference defaults';
if isfield(inputs,'params') && isfield(inputs,'time_s') && isfield(inputs,'irradiance_W_m2') && isfield(inputs,'cellTemperature_C')
    old=inputs; checkStruct(old.params,'arrayPower params');
    names=fieldnames(p);
    for k=1:numel(names)
        if isfield(old.params,names{k}), p.(names{k})=old.params.(names{k}); end
    end
    if isfield(old.params,'selectedArrayArea_m2'), p.area_m2=old.params.selectedArrayArea_m2; end
    inputs=struct('time_s',old.time_s,'irradiance_W_m2',old.irradiance_W_m2,'cellTemperature_C',old.cellTemperature_C);
    if isfield(old,'thermal') && isfield(old.thermal,'temperatureRate_K_s')
        inputs.temperatureRate_K_s=old.thermal.temperatureRate_K_s;
    end
    if isfield(old,'backTemperature_C'), inputs.backTemperature_C=old.backTemperature_C; end
    if isfield(old,'thermal') && isfield(old.thermal,'backTemperatureRate_K_s')
        inputs.backTemperatureRate_K_s=old.thermal.backTemperatureRate_K_s;
    end
    if isfield(old,'environment')
        for boundaryName={'frontEarthViewFactor','backEarthViewFactor','frontAbsorbedEarthFlux_W_m2','backAbsorbedEarthFlux_W_m2'}
            inputs.(boundaryName{1})=old.environment.(boundaryName{1});
        end
    end
    source='arrayPower settings, projected irradiance and supplied front/back temperature history';
end
names=fieldnames(overrides);
for k=1:numel(names)
    if ~isfield(p,names{k}), error('SolarArrayThermalModel:parameter','Unknown parameter: %s',names{k}); end
    p.(names{k})=overrides.(names{k});
end
validateSettings(p);
validateattributes(p.thermalNodeCount,{'numeric'},{'scalar','integer','>=',1,'<=',2});
allowed={'time_s','timeEdges_s','irradiance_W_m2','normalIrradiance_W_m2','incidenceAngle_deg', ...
    'illumination','cellTemperature_C','backTemperature_C','backTemperatureRate_K_s','temperatureRate_K_s','electricalPower_W','operatingPointFraction', ...
    'heaterPower_W','mountTemperature_C','structureTemperature_C','absorbedEarthFlux_W_m2', ...
    'frontEarthViewFactor','backEarthViewFactor','frontAbsorbedEarthFlux_W_m2','backAbsorbedEarthFlux_W_m2'};
names=fieldnames(inputs);
for k=1:numel(names)
    if ~ismember(names{k},allowed), error('SolarArrayThermalModel:input','Unknown input: %s',names{k}); end
end
if any(isfield(inputs,{'backTemperature_C','backTemperatureRate_K_s'})) && ~isfield(inputs,'cellTemperature_C')
    error('SolarArrayThermalModel:temperature','Supplied rear temperature/derivative requires a supplied front temperature.');
end
if isfield(inputs,'backTemperatureRate_K_s') && ~isfield(inputs,'backTemperature_C')
    error('SolarArrayThermalModel:temperature','A supplied rear derivative requires a supplied rear temperature.');
end
if p.thermalNodeCount==2 && xor(isfield(inputs,'temperatureRate_K_s'),isfield(inputs,'backTemperatureRate_K_s'))
    error('SolarArrayThermalModel:temperature','Two-node derivative postprocessing requires both front and rear derivatives.');
end
[t,edges,n]=timeBasis(inputs);
if isfield(inputs,'temperatureRate_K_s') && ~isfield(inputs,'cellTemperature_C')
    error('SolarArrayThermalModel:temperature','A supplied temperature derivative requires supplied temperature.');
end
if isfield(inputs,'irradiance_W_m2')
    if any(isfield(inputs,{'normalIrradiance_W_m2','incidenceAngle_deg','illumination'}))
        error('SolarArrayThermalModel:flux','Do not apply incidence/eclipse to already-projected irradiance.');
    end
    S=signal(inputs,'irradiance_W_m2',0,n);
else
    normal=signal(inputs,'normalIrradiance_W_m2',p.defaultNormalIrradiance_W_m2,n);
    angle=signal(inputs,'incidenceAngle_deg',p.defaultIncidenceAngle_deg,n);
    illumination=signal(inputs,'illumination',p.defaultIllumination,n);
    if any(normal<0 | angle<0 | angle>180 | illumination<0 | illumination>1)
        error('SolarArrayThermalModel:flux','Require nonnegative flux, angle in [0,180] and illumination in [0,1].');
    end
    S=normal.*max(cosd(angle),0).*illumination;
end
heater=signal(inputs,'heaterPower_W',p.defaultHeaterPower_W,n);
mountT=signal(inputs,'mountTemperature_C',p.mountTemperature_C,n);
structureT=signal(inputs,'structureTemperature_C',p.structureTemperature_C,n);
earth=signal(inputs,'absorbedEarthFlux_W_m2',p.absorbedEarthFlux_W_m2,n);
boundaries=struct();
for boundaryName={'frontEarthViewFactor','backEarthViewFactor','frontAbsorbedEarthFlux_W_m2','backAbsorbedEarthFlux_W_m2'}
    name=boundaryName{1}; boundaries.(name)=signal(inputs,name,p.(name),n);
    if any(boundaries.(name)<0), error('SolarArrayThermalModel:boundary','Earth fluxes and view factors must be nonnegative.'); end
end
if any(boundaries.frontEarthViewFactor+p.frontStructureViewFactor>1+1e-12) || ...
        any(boundaries.backEarthViewFactor+p.backStructureViewFactor>1+1e-12)
    error('SolarArrayThermalModel:view','Earth and structure view fractions cannot exceed one per face.');
end
if any(S<0 | heater<0 | earth<0 | mountT<=-273.15 | structureT<=-273.15)
    error('SolarArrayThermalModel:physics','Fluxes/heater must be nonnegative; sink temperatures must exceed absolute zero.');
end
actualExport=isfield(inputs,'electricalPower_W');
if actualExport && isfield(inputs,'operatingPointFraction')
    error('SolarArrayThermalModel:export','Supply actual export OR an operating-point fraction.');
end
if actualExport, electrical=signal(inputs,'electricalPower_W',0,n); else, electrical=[]; end
fraction=signal(inputs,'operatingPointFraction',p.defaultOperatingPointFraction,n);
if any(fraction<0 | fraction>1), error('SolarArrayThermalModel:export','Operating-point fraction must be in [0,1].'); end
absorbed=p.solarAbsorptivity*S*p.selfShadowRetained*p.area_m2;
if actualExport && any(electrical<0 | electrical>absorbed+p.balanceTolerance_W)
    error('SolarArrayThermalModel:export','Terminal export must lie between zero and absorbed solar power.');
end

%% Solve or postprocess the panel temperature
if ~isempty(edges) && ~isfield(inputs,'cellTemperature_C')
    [T_C,Tedges,values,backEdges]=solveIntervals(edges,S,electrical,fraction,heater,mountT,structureT,earth,p,boundaries);
    storage=p.area_m2*heatCapacity(p).*diff((Tedges+backEdges)/2)./diff(edges);
    temperatureMode='solved transient history using actual export or selected operating fraction';
    storageSource='total capacitance times mean front/rear boundary-temperature change / interval duration';
else
    if n>1 && ~isfield(inputs,'cellTemperature_C')
        error('SolarArrayThermalModel:temperature','Use timeEdges_s to solve a transient, or supply temperature for time_s samples.');
    end
    T_C=signal(inputs,'cellTemperature_C',p.defaultCellTemperature_C,n);
    if any(T_C<=-273.15), error('SolarArrayThermalModel:temperature','Panel temperature must exceed absolute zero.'); end
    back_C=signal(inputs,'backTemperature_C',T_C,n);
    if p.thermalNodeCount==1 && any(abs(back_C-T_C)>1e-9), error('SolarArrayThermalModel:temperature','One thermal node requires equal front/back temperatures.'); end
    if any(back_C<=-273.15), error('SolarArrayThermalModel:temperature','Back temperature must exceed absolute zero.'); end
    q=p; q.mountTemperature_C=mountT; q.structureTemperature_C=structureT; q.absorbedEarthFlux_W_m2=earth; q=applyBoundaries(q,boundaries,1:n);
    if actualExport
        h=finalPowerEPS.solarArrayHeatBalance(S,T_C+273.15,q,electrical/p.area_m2,heater/p.area_m2,back_C+273.15);
    else
        h=finalPowerEPS.solarArrayHeatBalance(S,T_C+273.15,q,[],0,back_C+273.15);
        h=finalPowerEPS.solarArrayHeatBalance(S,T_C+273.15,q,h.potentialElectrical_W_m2.*fraction,heater/p.area_m2,back_C+273.15);
    end
    values=heatVector(h,T_C+273.15,back_C+273.15); Tedges=[]; backEdges=[];
    if isfield(inputs,'temperatureRate_K_s')
        frontRate=signal(inputs,'temperatureRate_K_s',0,n);
        backRate=signal(inputs,'backTemperatureRate_K_s',frontRate,n);
        storage=p.area_m2*heatCapacity(p)*(frontRate+backRate)/2;
        storageSource='supplied temperature derivative; residual independently checks consistency';
    else
        storage=p.area_m2*h.storage_W_m2;
        storageSource='inferred from heat balance; no independent temperature derivative supplied';
    end
    temperatureMode='supplied temperature history or labeled reference point';
end

%% Heat rejected, heat deposited and thermal storage
A=p.area_m2;
result.time_s=t; result.timeEdges_s=edges;
result.QradiationToSpace_W=A*values(5,:);
result.QradiationToEarth_W=A*values(17,:);
result.environmentBoundaries=boundaries;
result.QradiationToStructure_W=A*values(6,:);
result.QconductionToStructure_W=A*values(7,:);
result.QtoStructure_W=result.QradiationToStructure_W+result.QconductionToStructure_W;
result.Qout_W=A*values(8,:); % PRIMARY: net total thermal power leaving panel
result.QheatInput_W=A*values(9,:);
result.Qstorage_W=storage;
result.QabsorbedSolar_W=A*values(3,:);
result.QabsorbedEarth_W=A*values(4,:);
result.Qheater_W=heater;
result.electricalPower_W=A*values(1,:);
result.potentialElectricalPower_W=A*values(2,:);
result.extractionExceedsAvailable=result.electricalPower_W>result.potentialElectricalPower_W+p.balanceTolerance_W;
result.cellTemperature_C=T_C; result.temperatureEdges_C=Tedges;
result.backTemperature_C=values(11,:)-273.15; result.backTemperatureEdges_C=backEdges;
result.meanPanelTemperature_C=(T_C+result.backTemperature_C)/2;
result.QfrontRadiation_W=A*values(12,:); result.QbackRadiation_W=A*values(13,:);
result.QfrontToBack_W=A*values(14,:);
result.meanTemperatureRate_K_s=storage/(A*heatCapacity(p));
result.temperatureRate_K_s=values(15,:); result.backTemperatureRate_K_s=values(16,:);
result.thermalCapacitance_J_K=A*heatCapacity(p);
result.balanceResidual_W=result.QheatInput_W-result.Qout_W-result.Qstorage_W;
result.Qradiation_W=result.QradiationToSpace_W+result.QradiationToEarth_W+result.QradiationToStructure_W;
result.QnetWithoutConduction_W=result.QheatInput_W-result.Qradiation_W;
result.QtotalHeat_W=result.QheatInput_W; % retained name: DEPOSITED heat, not rejected heat
result.Qgen_W=result.QheatInput_W; % legacy alias: total deposited heat, no longer wiring-only
result.QsolarHeat_W=result.QabsorbedSolar_W-result.electricalPower_W;
result.Qgen_W_m2=result.Qgen_W/A;
result.QoutAtIntervalStart_W=[]; result.QoutAtIntervalEnd_W=[];
if ~isempty(Tedges)
    q=p; q.mountTemperature_C=mountT; q.structureTemperature_C=structureT; q.absorbedEarthFlux_W_m2=earth; q=applyBoundaries(q,boundaries,1:n);
    left=boundaryHeat(Tedges(1:end-1)+273.15,backEdges(1:end-1)+273.15);
    right=boundaryHeat(Tedges(2:end)+273.15,backEdges(2:end)+273.15);
    result.intervalStart=boundaryPowers(left,A);
    result.intervalEnd=boundaryPowers(right,A);
    result.QoutAtIntervalStart_W=result.intervalStart.Qout_W;
    result.QoutAtIntervalEnd_W=result.intervalEnd.Qout_W;
end
result.irradiance_W_m2=S;
result.temperatureExtrapolated=(S>0) & (T_C<p.powerTempCoeffRange_C(1) | T_C>p.powerTempCoeffRange_C(2));
result.productTemperatureExceeded=T_C<p.productOperatingRange_C(1) | T_C>p.productOperatingRange_C(2) | ...
    result.backTemperature_C<p.productOperatingRange_C(1) | result.backTemperature_C>p.productOperatingRange_C(2);
if ~isempty(Tedges)
    result.productTemperatureExceeded=result.productTemperatureExceeded | ...
        Tedges(1:end-1)<p.productOperatingRange_C(1) | Tedges(1:end-1)>p.productOperatingRange_C(2) | ...
        Tedges(2:end)<p.productOperatingRange_C(1) | Tedges(2:end)>p.productOperatingRange_C(2) | ...
        backEdges(1:end-1)<p.productOperatingRange_C(1) | backEdges(1:end-1)>p.productOperatingRange_C(2) | ...
        backEdges(2:end)<p.productOperatingRange_C(1) | backEdges(2:end)>p.productOperatingRange_C(2);
end
result.summary=heatSummary(result.Qout_W,t,edges);
result.outflowSummary=result.summary;
result.spaceRadiationSummary=heatSummary(result.QradiationToSpace_W,t,edges);
result.structureHeatSummary=heatSummary(result.QtoStructure_W,t,edges);
result.totalHeatSummary=heatSummary(result.QheatInput_W,t,edges);
result.storageSummary=heatSummary(result.Qstorage_W,t,edges);
result.absorbedSolarSummary=heatSummary(result.QabsorbedSolar_W,t,edges);
if ~isempty(Tedges)
    % Mean/energy use interval integrals; extrema include instantaneous interval boundaries.
    result.outflowSummary=includeBounds(result.outflowSummary,result.intervalStart.Qout_W,result.intervalEnd.Qout_W);
    result.summary=result.outflowSummary;
    result.spaceRadiationSummary=includeBounds(result.spaceRadiationSummary,result.intervalStart.QradiationToSpace_W,result.intervalEnd.QradiationToSpace_W);
    result.structureHeatSummary=includeBounds(result.structureHeatSummary,result.intervalStart.QtoStructure_W,result.intervalEnd.QtoStructure_W);
    result.totalHeatSummary=includeBounds(result.totalHeatSummary,result.intervalStart.QheatInput_W,result.intervalEnd.QheatInput_W);
    result.storageSummary=includeBounds(result.storageSummary,result.intervalStart.Qstorage_W,result.intervalEnd.Qstorage_W);
end
screen=finalPowerEPS.solarArrayHeatBalance(S,T_C+273.15,q,[],0,result.backTemperature_C+273.15);
result.cellEfficiency=screen.cellEfficiency;
result.lumpedBiot=(screen.radiationCoefficient_W_m2K+p.mountConductancePerArea_W_m2K)* ...
    p.panelThickness_m/p.panelConductivity_W_mK;
result.maxLumpedBiot=max(result.lumpedBiot);
if ~isempty(Tedges)
    edgeCoefficient=max([left.radiationCoefficient_W_m2K,right.radiationCoefficient_W_m2K]);
    result.maxLumpedBiot=max(result.maxLumpedBiot,(edgeCoefficient+p.mountConductancePerArea_W_m2K)* ...
        p.panelThickness_m/p.panelConductivity_W_mK);
end
result.params=p; result.inputSource=source;
result.temperatureMode=temperatureMode; result.storageSource=storageSource;
result.structureExchangeEnabled=p.mountConductancePerArea_W_m2K>0 || ...
    p.frontEmissivity*p.frontStructureViewFactor+p.backEmissivity*p.backStructureViewFactor>0;
result.limitations=['Two half-plate thermal nodes resolve a front/rear temperature difference. Rear emissivity remains a proxy. Earth heating follows supplied orbital boundaries; unknown mounting/structure paths remain disabled. ' ...
    'A disabled-path zero is not evidence of zero physical transfer. Structure temperatures are prescribed; radiation assumes ' ...
    'large isothermal black surroundings, not an enclosure solution. Two nodes omit in-plane gradients and individual layers/hot spots. ' ...
    'Qgen/QtotalHeat are legacy deposited-heat aliases; Qout is heat leaving. Supplied temperatures are postprocessed.'];
if p.buildSummary, result.summaryTable=finalPowerEPS.solarHeatRateTable(result); else, result.summaryTable=[]; end
if any(~isfinite([result.Qout_W result.QheatInput_W result.Qstorage_W]))
    error('SolarArrayThermalModel:numerics','Nonfinite thermal result.');
end
if p.showSummary
    fprintf('\nSOLAR ARRAY THERMAL SUMMARY\n');
    if result.summary.duration_s==0
        fprintf('Solar array: supplied/reference operating point at %.4g C; no orbital temperature history solved.\n',T_C(1));
    end
    disp(result.summaryTable);
    fprintf(['Table key [W]:\n' ...
        '  Qout                Total net heat leaving through enabled radiation/conduction paths.\n' ...
        '  QradiationToSpace   Net radiation exchanged with space.\n' ...
        '  QtoStructure        Structure radiation + mount conduction; positive heats the spacecraft.\n' ...
        '  QheatInput          Absorbed solar/Earth - terminal electrical export + heater.\n' ...
        '  Qstorage            Rate of thermal energy storage; positive warming, negative cooling.\n\n' ...
        'Notes:\n' ...
        '  QheatInput = Qout + Qstorage. Positive transfer is outward; negative structure transfer heats the panel.\n' ...
        '  Means are time-weighted when a history is supplied. Space rejection is not an internal spacecraft heat load.\n']);
    if ~result.structureExchangeEnabled
        fprintf('  Structure paths are DISABLED pending conductance/view factors and sink temperatures; zero is not a verified physical load.\n');
    end
    fprintf('\n');
end
if p.warnAssumptions && any(result.extractionExceedsAvailable)
    warning('SolarArrayThermalModel:export','Supplied export exceeds modeled available power; inspect extractionExceedsAvailable.');
end
if p.warnAssumptions && any(abs(result.balanceResidual_W)>p.balanceTolerance_W)
    warning('SolarArrayThermalModel:balance','Supplied temperature derivative is inconsistent; inspect balanceResidual_W.');
end
if p.warnAssumptions && any(result.temperatureExtrapolated | result.productTemperatureExceeded)
    warning('SolarArrayThermalModel:temperature','Some temperatures leave the proxy coefficient/product range; inspect flags.');
end
if p.warnAssumptions && result.maxLumpedBiot>p.lumpedBiotWarningLimit
    warning('SolarArrayThermalModel:lumped','Boundary heat transfer exceeds the conservative lumped Biot screening limit; check a layered model.');
end
    function h=boundaryHeat(value,back)
        if actualExport
            h=finalPowerEPS.solarArrayHeatBalance(S,value,q,electrical/A,heater/A,back);
        else
            h=finalPowerEPS.solarArrayHeatBalance(S,value,q,[],0,back);
            h=finalPowerEPS.solarArrayHeatBalance(S,value,q,h.potentialElectrical_W_m2.*fraction,heater/A,back);
        end
    end
end

%% Transient interval integrator
function [meanT_C,edgeT_C,averages,backEdges]=solveIntervals(edges,S,P,fraction,heater,mountT,structureT,earth,p,boundaries)
% RK4 integrates temperature and all heat channels with identical weights.
% Signals are constant interval values. Thermal-stiffness-limited substeps ensure stable cooling.
n=numel(S); averages=zeros(17,n); edgeT_C=zeros(1,n+1); edgeT_C(1)=p.initialTemperature_C;
backEdges=edgeT_C;
if p.thermalNodeCount==2, backEdges(1)=p.initialBackTemperature_C; end
T=p.initialTemperature_C+273.15;
if p.thermalNodeCount==2, T=[T;p.initialBackTemperature_C+273.15]; end
C=heatCapacity(p); sigma=5.670374419e-8;
for k=1:n
    q=p; q.mountTemperature_C=mountT(k); q.structureTemperature_C=structureT(k); q.absorbedEarthFlux_W_m2=earth(k); q=applyBoundaries(q,boundaries,k);
    duration=edges(k+1)-edges(k); total=zeros(17,1);
    if isempty(P), extracted=[]; else, extracted=P(k)/p.area_m2; end
    qin=p.solarAbsorptivity*S(k)*p.selfShadowRetained+earth(k)+boundaries.frontAbsorbedEarthFlux_W_m2(k)+boundaries.backAbsorbedEarthFlux_W_m2(k)+heater(k)/p.area_m2;
    radiation=sigma*(p.frontEmissivity+p.backEmissivity);
    sinks=max([p.spaceTemperature_K,mountT(k)+273.15,structureT(k)+273.15]);
    if radiation>0, bound=max(max(T),sinks+(max(0,qin)/radiation)^0.25);
    elseif p.mountConductancePerArea_W_m2K>0, bound=max(max(T),sinks+max(0,qin)/p.mountConductancePerArea_W_m2K);
    else, bound=max(T)+max(0,qin)*duration/C; end
    slope=4*radiation*bound^3+p.mountConductancePerArea_W_m2K;
    if p.thermalNodeCount==2
        slope=2*slope+8*p.panelConductivity_W_mK/p.panelThickness_m;
    end
    step=p.thermalMaxStep_s;
    if slope>0
        factor=0.2; if p.thermalNodeCount==2, factor=0.8; end
        step=min(step,factor*C/slope);
    end
    pieces=max(1,ceil(duration/step));
    if pieces>100000, error('SolarArrayThermalModel:stiff','Too many thermal substeps; review heat capacity/boundaries.'); end
    step=duration/pieces;
    for j=1:pieces
        [a,va]=evaluate(T); [b,vb]=evaluate(T+step*a/2);
        [c,vc]=evaluate(T+step*b/2); [d,vd]=evaluate(T+step*c);
        T=T+step*(a+2*b+2*c+d)/6;
        if any(~isfinite(T) | T<=0), error('SolarArrayThermalModel:temperature','Unphysical transient temperature.'); end
        total=total+step*(va+2*vb+2*vc+vd)/6;
    end
    averages(:,k)=total/duration; edgeT_C(k+1)=T(1)-273.15; backEdges(k+1)=T(end)-273.15;
end
meanT_C=averages(10,:)-273.15;
    function [rate,v]=evaluate(value)
        if any(value<=0), error('SolarArrayThermalModel:temperature','Nonpositive RK stage temperature.'); end
        h=finalPowerEPS.solarArrayHeatBalance(S(k),value(1),q,extracted,heater(k)/p.area_m2,value(end),false);
        if isempty(extracted) && fraction(k)~=1
            h=finalPowerEPS.solarArrayHeatBalance(S(k),value(1),q,h.potentialElectrical_W_m2*fraction(k),heater(k)/p.area_m2,value(end),false);
        end
        rate=h.temperatureRate_K_s;
        if p.thermalNodeCount==2, rate=[rate;h.backTemperatureRate_K_s]; end
        v=heatVector(h,value(1),value(end));
    end
end

function C=heatCapacity(p)
C=p.panelDensity_kg_m3*p.panelThickness_m*p.panelSpecificHeat_J_kgK+p.additionalArealHeatCapacity_J_m2K;
end
function validateSettings(p)
fractions={'cellEfficiency','cellCoverage','selfShadowRetained','mismatchRetained','opticalRetained', ...
    'electricalDegradationRetained','defaultOperatingPointFraction','solarAbsorptivity','frontEmissivity', ...
    'backEmissivity','frontStructureViewFactor','backStructureViewFactor','defaultIllumination'};
for k=1:numel(fractions), validateattributes(p.(fractions{k}),{'numeric'},{'real','finite','scalar','>=',0,'<=',1}); end
positive={'area_m2','panelThickness_m','panelDensity_kg_m3','panelSpecificHeat_J_kgK', ...
    'panelConductivity_W_mK','lumpedBiotWarningLimit','thermalMaxStep_s'};
for k=1:numel(positive), validateattributes(p.(positive{k}),{'numeric'},{'real','finite','scalar','positive'}); end
nonnegative={'spaceTemperature_K','mountConductancePerArea_W_m2K','absorbedEarthFlux_W_m2', ...
    'additionalArealHeatCapacity_J_m2K','defaultNormalIrradiance_W_m2','defaultHeaterPower_W', ...
    'energyBalanceTolerance_W_m2','balanceTolerance_W'};
for k=1:numel(nonnegative), validateattributes(p.(nonnegative{k}),{'numeric'},{'real','finite','scalar','nonnegative'}); end
temperatures={'initialBackTemperature_C','initialTemperature_C','defaultCellTemperature_C','referenceTemperature_C','structureTemperature_C','mountTemperature_C'};
for k=1:numel(temperatures), validateattributes(p.(temperatures{k}),{'numeric'},{'real','finite','scalar','>',-273.15}); end
validateattributes(p.powerTempCoeff_perC,{'numeric'},{'real','finite','scalar'});
validateattributes(p.defaultIncidenceAngle_deg,{'numeric'},{'real','finite','scalar','>=',0,'<=',180});
for name={'powerTempCoeffRange_C','productOperatingRange_C'}
    x=p.(name{1}); validateattributes(x,{'numeric'},{'real','finite','vector','numel',2,'>',-273.15});
    if x(2)<=x(1), error('SolarArrayThermalModel:parameter','Temperature ranges must increase.'); end
end
for name={'showSummary','warnAssumptions','buildSummary'}
    x=p.(name{1}); if ~(isnumeric(x)||islogical(x)) || ~isscalar(x) || ~isreal(x) || ~(x==0||x==1)
        error('SolarArrayThermalModel:parameter','Display flags must be true/false.');
    end
end
end
function checkStruct(s,name)
if ~isstruct(s)||~isscalar(s), error('SolarArrayThermalModel:input','%s must be a scalar struct.',name); end
end
function [t,edges,n]=timeBasis(s)
edges=[];
if isfield(s,'time_s') && isfield(s,'timeEdges_s'), error('SolarArrayThermalModel:time','Supply samples OR interval edges.'); end
if isfield(s,'timeEdges_s')
    edges=s.timeEdges_s; validateattributes(edges,{'numeric'},{'real','finite','vector','nonempty'}); edges=double(edges(:).');
    if numel(edges)<2 || any(diff(edges)<=0), error('SolarArrayThermalModel:time','Need two or more increasing interval edges.'); end
    t=edges(1:end-1)+diff(edges)/2;
else
    if isfield(s,'time_s'), t=s.time_s; else, t=0; end
    validateattributes(t,{'numeric'},{'real','finite','vector','nonempty'}); t=double(t(:).');
    if any(diff(t)<=0), error('SolarArrayThermalModel:time','Sample times must increase.'); end
end
n=numel(t);
end
function x=signal(s,name,default,n)
if isfield(s,name), x=s.(name); else, x=default; end
validateattributes(x,{'numeric'},{'real','finite','vector','nonempty'}); x=double(x(:).');
if isscalar(x), x=repmat(x,1,n); end
if numel(x)~=n, error('SolarArrayThermalModel:length','%s must be scalar or match the time basis.',name); end
end
function s=heatSummary(q,t,edges)
s.min_W=min(q); s.max_W=max(q);
if ~isempty(edges), s.duration_s=edges(end)-edges(1); s.heatEnergy_J=sum(q.*diff(edges));
elseif numel(t)>1, s.duration_s=t(end)-t(1); s.heatEnergy_J=trapz(t,q);
else, s.duration_s=0; s.heatEnergy_J=NaN; end
if s.duration_s>0, s.mean_W=s.heatEnergy_J/s.duration_s; else, s.mean_W=q(1); end
s.heatEnergy_Wh=s.heatEnergy_J/3600;
end

function s=includeBounds(s,left,right)
s.min_W=min([s.min_W,left,right]); s.max_W=max([s.max_W,left,right]);
end
function b=boundaryPowers(h,A)
b.Qout_W=A*h.heatOut_W_m2;
b.QradiationToSpace_W=A*h.radiationToSpace_W_m2;
b.QradiationToEarth_W=A*h.radiationToEarth_W_m2;
b.QtoStructure_W=A*h.heatToStructure_W_m2;
b.QheatInput_W=A*h.heatInput_W_m2;
b.Qstorage_W=A*h.storage_W_m2;
end

%% References and boundary definitions
% Cengel et al., Fundamentals of Thermal-Fluid Sciences, 5th ed., uploaded Ch.18 Eq.18-1 (storage)
% and Ch.21 Eq.21-53 (radiation to large surroundings). Mount conduction: Q=G*(Tpanel-Tmount).
% Two half-plate temperatures omit in-plane hot spots and individual adhesive/contact layers.
% Product data and optical/material proxies are documented in arrayPower.m.
% No wiring resistance or 4% resistive-loss assumption is needed for this total heat balance.

function p=applyBoundaries(p,boundaries,index)
for boundaryName={'frontEarthViewFactor','backEarthViewFactor','frontAbsorbedEarthFlux_W_m2','backAbsorbedEarthFlux_W_m2'}
    name=boundaryName{1}; p.(name)=boundaries.(name)(index);
end
end

function v=heatVector(h,T,back)
v=[h.electrical_W_m2;h.potentialElectrical_W_m2;h.absorbedSolar_W_m2;h.absorbedEarth_W_m2; ...
    h.radiationToSpace_W_m2;h.radiationToStructure_W_m2;h.conductionToStructure_W_m2; ...
    h.heatOut_W_m2;h.heatInput_W_m2;T;back;h.frontRadiation_W_m2;h.backRadiation_W_m2; ...
    h.frontToBack_W_m2;h.temperatureRate_K_s;h.backTemperatureRate_K_s;h.radiationToEarth_W_m2];
end
