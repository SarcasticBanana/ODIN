function result = componentThermalOrbit(out, cfg, arrayResult, overrides)
% COMPONENTTHERMALORBIT Nominal/peak component heat histories on the propagated orbit.
% Couples panel temperature/extraction to routed power, then reports total panel heat rejection.
% Usage: thermalResult = finalPowerEPS.componentThermalOrbit(out, cfg, arrayResult); main_demo does this automatically.

% Shared defaults and assumptions: edit finalPowerEPS.modelSettings('component').
p=finalPowerEPS.modelSettings('component');

if nargin == 0, result.params = p; return; end
if nargin < 3, error('componentThermalOrbit:input','Supply out, cfg and arrayResult from the same propagation.'); end
if nargin < 4, overrides = struct(); end
checkStruct(overrides, 'overrides');
names = fieldnames(overrides);
for k = 1:numel(names)
    if ~isfield(p,names{k}), error('componentThermalOrbit:parameter','Unknown parameter: %s',names{k}); end
    p.(names{k}) = overrides.(names{k});
end
validateSettings(p);
checkStruct(out,'out'); checkStruct(cfg,'cfg'); checkStruct(arrayResult,'arrayResult');
if ~isfield(out,'t') || ~isfield(out,'T') || ~isfield(arrayResult,'time_s')
    error('componentThermalOrbit:orbit','Require propagate output out.t/out.T and arrayPower time history.');
end
ts = double(arrayResult.time_s(:).');
orbitTime = double(out.t(:).');
validateattributes(ts,{'numeric'},{'real','finite','vector','nonempty'});
validateattributes(out.T,{'numeric'},{'real','finite','scalar','positive'});
if numel(ts)<2 || any(diff(ts)<=0) || ~isequal(size(ts),size(orbitTime)) || ...
        any(abs(ts-orbitTime)>1e-6)
    error('componentThermalOrbit:time','Array and propagated orbit must have matching increasing times.');
end
period = out.T;
if p.downlinkDuration_s > period/p.downlinkWindowsPerOrbit
    error('componentThermalOrbit:contacts','Downlink windows overlap; reduce duration or number of windows.');
end
if isfield(out,'illum')
    illum = double(out.illum(:).');
elseif isfield(cfg,'use') && isfield(cfg.use,'eclipse') && ~cfg.use.eclipse
    illum = ones(size(ts));
else
    error('componentThermalOrbit:illumination','Missing propagated illumination history; rerun propagate with eclipses enabled.');
end
if numel(illum)~=numel(ts) || any(~isfinite(illum) | illum<0 | illum>1)
    error('componentThermalOrbit:illumination','out.illum must match orbit samples and lie in [0,1].');
end

%% Shared interval grid and operation schedule
% Native orbit samples preserve the supplied eclipse boundaries; source flux is held until its next sample.
% Contact/camera boundaries are inserted explicitly, so the selected contact duration is not rounded.
start = ts(1); finish = ts(end);
edges = unique([ts, start:p.dispatchMaxStep_s:finish, finish, ...
    scheduleEdges(start,finish,period,p)]);
dt = diff(edges); mid = edges(1:end-1)+dt/2;
phase = mod((mid-start)/period,1);
contact = mod(phase-p.firstDownlinkPhase,1/p.downlinkWindowsPerOrbit) < p.downlinkDuration_s/period;
interpolation = 'linear';
if isfield(arrayResult,'params') && isfield(arrayResult.params,'irradianceInterpolation')
    interpolation = arrayResult.params.irradianceInterpolation;
end
sunFraction = interp1(ts,illum,mid,interpolation);
camera = mod(phase-p.cameraStartPhase,1) < p.cameraDutyFraction;
if p.cameraSunlitOnly, camera = camera & sunFraction>0; end
orbitIndex = floor((mid-start)/period)+1;

%% Component settings and available solar power
% Inherit the user's selected area (including 0.175 m^2 if changed in arrayPower), rather than imposing a new area.
% arrayPower temperatures initialize the iteration; nominal/peak each solve their own temperature history.
quiet = struct('showSummary',false,'warnAssumptions',false);
solarChanges = merge(p.solarOverrides,quiet);
if p.usePlanningAgingAllowance && isfield(arrayResult,'params') && ...
        isfield(arrayResult.params,'agingRetained') && ~isfield(p.solarOverrides,'electricalDegradationRetained')
    solarChanges.electricalDegradationRetained = arrayResult.params.agingRetained;
end
source = finalPowerEPS.SolarArrayThermalModel(arrayResult,solarChanges);
sp = source.params;
pmadChanges = p.pmadOverrides;
if isfield(arrayResult,'params')
    if isfield(arrayResult.params,'mpptEfficiency') && ~isfield(pmadChanges,'mpptEfficiency')
        pmadChanges.mpptEfficiency = arrayResult.params.mpptEfficiency;
    end
    if isfield(arrayResult.params,'converterEfficiency') && ~isfield(pmadChanges,'solarConverterEfficiency')
        pmadChanges.solarConverterEfficiency = arrayResult.params.converterEfficiency;
    end
end
pmadProbe = finalPowerEPS.PMADThermalModel(struct('solarToLoad_W',0,'batteryToLoad_W',0,'batteryCharge_W',0),merge(pmadChanges,quiet));
batteryProbe = finalPowerEPS.BatteryThermalModel(struct('packCurrent_A',0),merge(p.batteryOverrides,quiet));
mp = pmadProbe.params; bp = batteryProbe.params;
irradiance = interp1(ts,source.irradiance_W_m2,mid,interpolation);
temperatureGuess_C = interp1(ts,source.cellTemperature_C,mid,'linear');
solarInputs = struct('timeEdges_s',edges,'irradiance_W_m2',irradiance);
if isfield(arrayResult,'environment')
    for boundaryName={'frontEarthViewFactor','backEarthViewFactor','frontAbsorbedEarthFlux_W_m2','backAbsorbedEarthFlux_W_m2'}
        name=boundaryName{1};
        solarInputs.(name)=interp1(ts,arrayResult.environment.(name),mid,interpolation);
    end
end
sp.buildSummary=false; % summary tables are produced only after convergence

harnessEta = p.externalHarnessRetained;
if isnan(harnessEta)
    harnessEta = 1;
    if isfield(arrayResult,'params') && isfield(arrayResult.params,'wiringDiodeRetained')
        harnessEta = arrayResult.params.wiringDiodeRetained;
    end
end
validateattributes(harnessEta,{'numeric'},{'real','finite','scalar','>',0,'<=',1});
vocCell = p.cellOpenCircuitVoltage_V;
if isnan(vocCell), vocCell = bp.cellNominalVoltage_V; end
Tbattery = p.batteryCellTemperature_C;
if isnan(Tbattery), Tbattery = bp.defaultCellTemperature_C; end
if vocCell < bp.cellDischargeCutoffVoltage_V || vocCell > bp.cellChargeVoltage_V
    error('componentThermalOrbit:OCV','Assumed cell OCV must lie within cell voltage limits.');
end
chargeLimit = p.packChargeCurrentLimit_A;
if isnan(chargeLimit), chargeLimit = bp.parallelCells*bp.cellStandardChargeCurrent_A; end
chargeLimit = min(chargeLimit,bp.parallelCells*bp.cellMaxChargeCurrent_A);

%% Nominal and peak thermal sections
% Each scenario starts at the same SOC and uses the same sun/shadow/contact history.
% Peak means the upper rated powers while equipment is active, not an asserted transient duration.
labels = {'nominal','peak'};
for c = 1:2
    if c == 1
        load_W = p.computerNominal_W+p.adcsNominal_W+p.otherNominal_W + ...
            p.cameraNominal_W.*camera+p.commsNominal_W.*contact+p.commsStandby_W.*~contact;
    else
        load_W = p.computerPeak_W+p.adcsPeak_W+p.otherPeak_W + ...
            p.cameraPeak_W.*camera+p.commsPeak_W.*contact+p.commsStandby_W.*~contact;
    end
    guess = temperatureGuess_C;
    converged = false;
    for iteration=1:p.solarCouplingMaxIterations
        availableKernel=finalPowerEPS.solarArrayElectricalPower(irradiance,guess+273.15,sp);
        available.potentialElectricalPower_W=sp.area_m2*availableKernel.potentialElectrical_W_m2;
        d = routePower(load_W+sp.defaultHeaterPower_W,available.potentialElectricalPower_W,dt,p,mp,bp,harnessEta,vocCell,Tbattery,chargeLimit);
        suppliedBus = d.solarToLoad_W+d.batteryToLoad_W;
        served = min(load_W,suppliedBus); % active loads take priority over an optional panel heater
        modeSolar = solarInputs; % no supplied temperature: solve the thermal history
        modeSolar.electricalPower_W = d.panelExport_W;
        modeSolar.heaterPower_W = max(0,suppliedBus-served);
        a = finalPowerEPS.SolarArrayThermalModel(modeSolar,sp);
        temperatureError_C = max(abs(a.cellTemperature_C-guess));
        extractionError_W = max(abs(a.potentialElectricalPower_W-available.potentialElectricalPower_W));
        if temperatureError_C<=p.solarCouplingTolerance_C && extractionError_W<=p.solarCouplingTolerance_W && ...
                all(d.panelExport_W<=a.potentialElectricalPower_W+p.balanceTolerance_W)
            converged = true; break
        end
        guess = a.cellTemperature_C;
    end
    if ~converged
        error('componentThermalOrbit:solarCoupling','Panel temperature and power routing did not converge; review solar coupling settings.');
    end
    a.summaryTable=finalPowerEPS.solarHeatRateTable(a);
    a.params.buildSummary=true;
    mIn = struct('timeEdges_s',edges,'solarToLoad_W',d.solarToLoad_W,'batteryToLoad_W',d.batteryToLoad_W, ...
        'batteryCharge_W',d.batteryCharge_W,'controlPower_W',d.controlPower_W, ...
        'shuntPower_W',d.shuntPower_W,'extraSwitchCurrent_A',d.extraSwitchCurrent_A);
    m = finalPowerEPS.PMADThermalModel(mIn,mp);
    bIn = struct('timeEdges_s',edges,'packCurrent_A',d.packCurrent_A,'cellVoltage_V',d.cellVoltage_V, ...
        'cellTemperature_C',Tbattery,'heaterPower_W',d.batteryHeater_W,'electronicsPower_W',d.batteryElectronics_W);
    b = finalPowerEPS.BatteryThermalModel(bIn,bp);
    freeEnergyRate = -d.packCurrent_A.*(bp.seriesCells*vocCell);
    internalEnergyRate = freeEnergyRate-b.Qreversible_W;
    residual = a.QabsorbedSolar_W+a.QabsorbedEarth_W-served-a.QheatInput_W-d.harnessHeat_W- ...
        m.Qgen_W-b.Qgen_W-internalEnergyRate;
    if any(abs(residual)>p.balanceTolerance_W) || ...
            any(abs(b.packTerminalPower_W-d.externalBatteryPower_W)>p.balanceTolerance_W) || ...
            any(abs(a.balanceResidual_W)>p.balanceTolerance_W)
        error('componentThermalOrbit:balance','Cross-component energy accounting did not close.');
    end
    mode.load_W = load_W;
    mode.servedLoad_W = served;
    mode.unservedLoad_W = max(0,load_W-served);
    mode.solar = a; mode.pmad = m; mode.battery = b;
    mode.solarCouplingIterations = iteration;
    mode.solarCouplingTemperatureError_C = temperatureError_C;
    mode.solarCouplingPowerError_W = extractionError_W;
    mode.batterySOC = d.SOC;
    mode.harnessHeat_W = d.harnessHeat_W;
    mode.unservedAuxiliary_W = d.unservedAuxiliary_W+max(0,sp.defaultHeaterPower_W-a.Qheater_W);
    mode.panelCurtailment_W = max(0,available.potentialElectricalPower_W-d.panelExport_W);
    mode.internalChemicalEnergyRate_W = internalEnergyRate;
    mode.energyResidual_W = residual;
    mode.unservedEnergy_Wh = sum(mode.unservedLoad_W.*dt)/3600;
    mode.internalChemicalEnergyChange_Wh = sum(internalEnergyRate.*dt)/3600;
    mode.perOrbit = orbitSummary(mode,dt,orbitIndex,sunFraction);
    result.(labels{c}) = mode;
end

%% Results and command-window key
result.time_s = mid;
result.timeEdges_s = edges;
result.intervalDuration_s = dt;
result.orbitPeriod_s = period;
result.orbitIndex = orbitIndex;
result.illumination = sunFraction;
result.irradianceInterpolation = interpolation;
result.downlinkActive = contact;
result.cameraActive = camera;
result.params = p;
result.assumedCellOCV_V = vocCell;
result.harnessRetained = harnessEta;
result.packChargeCurrentLimit_A = chargeLimit;
result.effectiveCapacity_Ah = bp.parallelCells*bp.cellNominalCapacity_Ah*p.capacityRetained;
result.limitations = ['One propagated spacecraft history, not a sum of constellation spacecraft. ' ...
    'Contacts are scheduled assumptions; actual ground-station visibility is not calculated. ' ...
    'Blank component loads/standby powers remain unbudgeted. Constant OCV and assumed resistance replace ' ...
    'unavailable SOC/temperature-dependent cell curves. Panel temperature is coupled to actual extraction; ' ...
    'battery and PMAD temperatures are not solved here. Structure heat paths need measured/estimated boundaries.'];
result.summaryTables = thermalSummaryTables(result,sp,bp,mp,Tbattery);
if p.showSummary
    printThermalReport(result);
end
if any(result.nominal.solar.temperatureExtrapolated | result.nominal.solar.productTemperatureExceeded) || ...
        any(result.peak.solar.temperatureExtrapolated | result.peak.solar.productTemperatureExceeded)
    warning('componentThermalOrbit:arrayTemperature', ...
        'A solved panel temperature leaves the proxy coefficient or operating range; inspect solar temperature flags and the retained hot initial condition.');
end
if result.nominal.solar.maxLumpedBiot>result.nominal.solar.params.lumpedBiotWarningLimit || ...
        result.peak.solar.maxLumpedBiot>result.peak.solar.params.lumpedBiotWarningLimit
    warning('componentThermalOrbit:lumpedModel','An enabled boundary path exceeds the conservative lumped-model screening limit; review a layered panel model.');
end
if any(result.nominal.unservedLoad_W>p.balanceTolerance_W) || any(result.peak.unservedLoad_W>p.balanceTolerance_W) || ...
        any(result.nominal.unservedAuxiliary_W>p.balanceTolerance_W) || any(result.peak.unservedAuxiliary_W>p.balanceTolerance_W)
    warning('componentThermalOrbit:unserved','Some loads cannot be supplied at the imposed SOC/current/voltage limits; see unservedLoad_W.');
end
end

%% Power routing and battery-current section
function d = routePower(load_W,panelAvailable,dt,p,mp,bp,harnessEta,vocCell,Tbattery,chargeLimit)
n = numel(dt); zero=zeros(1,n);
d.solarToLoad_W=zero; d.batteryToLoad_W=zero; d.batteryCharge_W=zero;
d.panelExport_W=zero; d.packCurrent_A=zero; d.cellVoltage_V=zero;
d.controlPower_W=zero; d.shuntPower_W=zero; d.extraSwitchCurrent_A=zero;
d.batteryHeater_W=zero; d.batteryElectronics_W=zero; d.unservedAuxiliary_W=zero;
d.externalBatteryPower_W=zero; d.harnessHeat_W=zero;
d.SOC=zeros(1,n+1); d.SOC(1)=p.initialSOC;
etaS=mp.mpptEfficiency*mp.solarConverterEfficiency;
etaC=mp.mpptEfficiency*mp.batteryChargeConverterEfficiency;
etaD=mp.batteryDischargeEfficiency;
RcellPack=bp.seriesCells/bp.parallelCells*bp.cellDCResistance_Ohm;
Rpack=RcellPack+bp.packInterconnectResistance_Ohm;
Voc=bp.seriesCells*vocCell;
capacity_As=bp.parallelCells*bp.cellNominalCapacity_Ah*p.capacityRetained*3600;
auxM=mp.defaultControlPower_W+mp.defaultShuntPower_W+mp.defaultExtraSwitchCurrent_A^2*mp.extraSwitchResistance_Ohm;
auxB=bp.defaultHeaterPower_W+bp.defaultElectronicsPower_W;
maxDischarge=bp.parallelCells*bp.cellMaxDischargeCurrent_A;
maxCharge=chargeLimit;
if Rpack>0
    maxDischarge=min(maxDischarge,Voc/(2*Rpack)); % stable power-transfer branch
    maxDischarge=min(maxDischarge,max(0,(Voc-bp.seriesCells*bp.cellDischargeCutoffVoltage_V)/Rpack));
    maxCharge=min(maxCharge,max(0,(bp.seriesCells*bp.cellChargeVoltage_V-Voc)/Rpack));
end
if Tbattery<bp.cellDischargeTemperatureRange_C(1) || Tbattery>bp.cellDischargeTemperatureRange_C(2), maxDischarge=0; end
if Tbattery<bp.cellChargeTemperatureRange_C(1) || Tbattery>bp.cellChargeTemperatureRange_C(2), maxCharge=0; end
for k=1:n
    PV=panelAvailable(k)*harnessEta;
    auxSolar=min(PV,auxM);
    solarLoad=min(load_W(k),max(0,PV-auxSolar)*etaS);
    deficit=max(0,load_W(k)-solarLoad);
    auxNeed=auxM-auxSolar;
    Imax=min(maxDischarge,max(0,(d.SOC(k)-p.minimumSOC)*capacity_As/dt(k)));
    Jmax=min(maxCharge,max(0,(p.maximumSOC-d.SOC(k))*capacity_As/dt(k)));
    if deficit>0 || auxNeed>0
        request=deficit/etaD+auxNeed+auxB;
        I=min(Imax,currentFromNetworkPower(request,Voc,Rpack));
        networkPower=I*(Voc-I*Rpack);
        auxBattery=min(auxB,networkPower);
        external=max(0,networkPower-auxBattery);
        auxBatteryToPMAD=min(auxNeed,external);
        batteryLoad=min(deficit,max(0,external-auxBatteryToPMAD)*etaD);
        charge=0;
        actualAuxM=auxSolar+auxBatteryToPMAD;
        pmadSolarInput=auxSolar+solarLoad/etaS;
    else
        chargeAvailable=max(0,PV-auxSolar-solarLoad/etaS)*etaC;
        request=auxB-chargeAvailable;
        I=currentFromNetworkPower(request,Voc,Rpack);
        I=min(Imax,max(-Jmax,I));
        networkPower=I*(Voc-I*Rpack);
        auxBattery=min(auxB,max(0,chargeAvailable+networkPower));
        charge=max(0,auxBattery-networkPower);
        batteryLoad=0;
        actualAuxM=auxM;
        pmadSolarInput=auxSolar+solarLoad/etaS+charge/etaC;
    end
    fM=1; if auxM>0, fM=actualAuxM/auxM; end
    fB=1; if auxB>0, fB=auxBattery/auxB; end
    d.solarToLoad_W(k)=solarLoad; d.batteryToLoad_W(k)=batteryLoad; d.batteryCharge_W(k)=charge;
    d.panelExport_W(k)=pmadSolarInput/harnessEta;
    d.harnessHeat_W(k)=d.panelExport_W(k)-pmadSolarInput;
    d.packCurrent_A(k)=I;
    d.cellVoltage_V(k)=(Voc-I*RcellPack)/bp.seriesCells;
    d.externalBatteryPower_W(k)=networkPower-auxBattery;
    d.controlPower_W(k)=mp.defaultControlPower_W*fM;
    d.shuntPower_W(k)=mp.defaultShuntPower_W*fM;
    d.extraSwitchCurrent_A(k)=mp.defaultExtraSwitchCurrent_A*sqrt(max(0,fM));
    d.batteryHeater_W(k)=bp.defaultHeaterPower_W*fB;
    d.batteryElectronics_W(k)=bp.defaultElectronicsPower_W*fB;
    d.unservedAuxiliary_W(k)=max(0,auxM-actualAuxM)+max(0,auxB-auxBattery);
    d.SOC(k+1)=d.SOC(k)-I*dt(k)/capacity_As;
end
if any(d.SOC<p.minimumSOC-1e-10 | d.SOC>p.maximumSOC+1e-10) || any(d.panelExport_W>panelAvailable+p.balanceTolerance_W)
    error('componentThermalOrbit:dispatch','SOC or available solar power bounds were violated.');
end
end

function I = currentFromNetworkPower(P,Voc,R)
% P=I*(Voc-I*R); positive discharges, negative charges. Clamp an infeasible demand to peak transfer before current limiting.
D=max(0,Voc^2-4*R*P);
I=2*P/(Voc+sqrt(D));
if P>0 && R>0 && P>Voc^2/(4*R), I=Voc/(2*R); end
end

%% Per-orbit summaries and input helpers
function rows = orbitSummary(mode,dt,index,illum)
for j=1:max(index)
    use=index==j; duration=sum(dt(use));
    rows(j).Orbit=j; rows(j).Duration_min=duration/60;
    rows(j).Umbra_min=sum(dt(use & illum==0))/60;
    rows(j).Illuminated_min=sum(dt(use & illum>0))/60;
    rows(j).MeanPMAD_W=sum(mode.pmad.Qgen_W(use).*dt(use))/duration;
    rows(j).MinPMAD_W=min(mode.pmad.Qgen_W(use)); rows(j).MaxPMAD_W=max(mode.pmad.Qgen_W(use));
    rows(j).MeanBattery_W=sum(mode.battery.Qgen_W(use).*dt(use))/duration;
    rows(j).MinBattery_W=min(mode.battery.Qgen_W(use)); rows(j).MaxBattery_W=max(mode.battery.Qgen_W(use));
    rows(j).MeanPanelHeating_W=sum(mode.solar.QtotalHeat_W(use).*dt(use))/duration;
    rows(j).MinPanelHeating_W=min(mode.solar.QtotalHeat_W(use)); rows(j).MaxPanelHeating_W=max(mode.solar.QtotalHeat_W(use));
    rows(j).MeanSolarOutput_W=sum(mode.solar.Qout_W(use).*dt(use))/duration;
    rows(j).MinSolarOutput_W=min([mode.solar.QoutAtIntervalStart_W(use),mode.solar.QoutAtIntervalEnd_W(use)]);
    rows(j).MaxSolarOutput_W=max([mode.solar.QoutAtIntervalStart_W(use),mode.solar.QoutAtIntervalEnd_W(use)]);
    rows(j).MeanSolarToSpace_W=sum(mode.solar.QradiationToSpace_W(use).*dt(use))/duration;
    rows(j).MeanSolarToStructure_W=sum(mode.solar.QtoStructure_W(use).*dt(use))/duration;
    rows(j).MeanPanelStorage_W=sum(mode.solar.Qstorage_W(use).*dt(use))/duration;
    rows(j).MinSOC_percent=100*min(mode.batterySOC([find(use),find(use,1,'last')+1]));
    rows(j).UnservedLoad_Wh=sum(mode.unservedLoad_W(use).*dt(use))/3600;
end
end
function edges=scheduleEdges(start,finish,period,p)
edges=[];
for orbit=-1:ceil((finish-start)/period)
    contacts=start+(orbit+p.firstDownlinkPhase+(0:p.downlinkWindowsPerOrbit-1)/p.downlinkWindowsPerOrbit)*period;
    camera=start+(orbit+p.cameraStartPhase)*period;
    edges=[edges,contacts,contacts+p.downlinkDuration_s,camera,camera+p.cameraDutyFraction*period,start+orbit*period]; %#ok<AGROW>
end
edges=edges(edges>start & edges<finish);
end
function validateSettings(p)
nonnegative={'computerNominal_W','computerPeak_W','adcsNominal_W','adcsPeak_W','cameraNominal_W', ...
    'cameraPeak_W','commsNominal_W','commsPeak_W','commsStandby_W','otherNominal_W','otherPeak_W','downlinkDuration_s','balanceTolerance_W'};
for k=1:numel(nonnegative), validateattributes(p.(nonnegative{k}),{'numeric'},{'real','finite','scalar','nonnegative'}); end
validateattributes(p.downlinkWindowsPerOrbit,{'numeric'},{'real','finite','scalar','integer','positive'});
validateattributes(p.dispatchMaxStep_s,{'numeric'},{'real','finite','scalar','positive'});
validateattributes(p.solarCouplingMaxIterations,{'numeric'},{'real','finite','scalar','integer','positive'});
validateattributes(p.solarCouplingTolerance_C,{'numeric'},{'real','finite','scalar','positive'});
validateattributes(p.solarCouplingTolerance_W,{'numeric'},{'real','finite','scalar','positive'});
unit={'cameraDutyFraction','initialSOC','minimumSOC','maximumSOC','capacityRetained'};
for k=1:numel(unit), validateattributes(p.(unit{k}),{'numeric'},{'real','finite','scalar','>=',0,'<=',1}); end
phase={'firstDownlinkPhase','cameraStartPhase'};
for k=1:numel(phase), validateattributes(p.(phase{k}),{'numeric'},{'real','finite','scalar','>=',0,'<',1}); end
if p.minimumSOC>p.initialSOC || p.initialSOC>p.maximumSOC || p.capacityRetained<=0 || p.minimumSOC>=p.maximumSOC
    error('componentThermalOrbit:SOC','Require minimumSOC<=initialSOC<=maximumSOC and positive capacity.');
end
if any([p.computerPeak_W p.adcsPeak_W p.cameraPeak_W p.commsPeak_W p.otherPeak_W] < ...
        [p.computerNominal_W p.adcsNominal_W p.cameraNominal_W p.commsNominal_W p.otherNominal_W])
    error('componentThermalOrbit:load','Peak powers must not be below nominal powers.');
end
optional={'cellOpenCircuitVoltage_V','batteryCellTemperature_C','packChargeCurrentLimit_A','externalHarnessRetained'};
for k=1:numel(optional)
    x=p.(optional{k}); validateattributes(x,{'numeric'},{'real','scalar'});
    if ~isnan(x) && ~isfinite(x), error('componentThermalOrbit:parameter','Optional settings must be finite or NaN.'); end
end
if ~isnan(p.packChargeCurrentLimit_A) && p.packChargeCurrentLimit_A<0, error('componentThermalOrbit:parameter','Charge current limit must be nonnegative.'); end
if ~isnan(p.batteryCellTemperature_C) && p.batteryCellTemperature_C<=-273.15, error('componentThermalOrbit:temperature','Battery temperature must exceed absolute zero.'); end
flags={'cameraSunlitOnly','usePlanningAgingAllowance','showSummary','showPerOrbitTables'};
for k=1:numel(flags)
    x=p.(flags{k});
    if ~(islogical(x)||isnumeric(x)) || ~isscalar(x) || ~isreal(x) || ~(x==0||x==1)
        error('componentThermalOrbit:parameter','%s must be true/false or 0/1.',flags{k});
    end
end
checkStruct(p.solarOverrides,'solarOverrides'); checkStruct(p.pmadOverrides,'pmadOverrides'); checkStruct(p.batteryOverrides,'batteryOverrides');
end
function checkStruct(x,name)
if ~isstruct(x)||~isscalar(x), error('componentThermalOrbit:input','%s must be a scalar struct.',name); end
end
function a=merge(a,b)
names=fieldnames(b); for k=1:numel(names), a.(names{k})=b.(names{k}); end
end
function tables=thermalSummaryTables(r,sp,bp,mp,Tbattery)
% Tables contain values only; definitions and modeling notes are printed separately.
p=r.params;
labels={'nominal','peak'};
keys={'PMAD_Qgen';'Battery_Qgen';'Solar_Qout';'Solar_Qspace'; ...
    'Solar_QtoStructure';'Solar_QheatInput';'Solar_Qstorage'};
Case=cell(14,1); Quantity=cell(14,1);
Mean_W=zeros(14,1); Min_W=Mean_W; Max_W=Mean_W;
LoadMin_W=zeros(2,1); LoadMax_W=LoadMin_W; MinimumSOC_percent=LoadMin_W;
Unserved_Wh=LoadMin_W; Min_C=LoadMin_W; Max_C=LoadMin_W; BackMin_C=LoadMin_W; BackMax_C=LoadMin_W; CouplingIterations=LoadMin_W;
for c=1:2
    v=r.(labels{c}); rows=(c-1)*7+(1:7);
    Case(rows)=repmat(labels(c),7,1); Quantity(rows)=keys;
    Mean_W(rows)=[v.pmad.summary.meanQgen_W;v.battery.summary.meanQgen_W; ...
        v.solar.outflowSummary.mean_W;v.solar.spaceRadiationSummary.mean_W; ...
        v.solar.structureHeatSummary.mean_W;v.solar.totalHeatSummary.mean_W;v.solar.storageSummary.mean_W];
    Min_W(rows)=[v.pmad.summary.minQgen_W;v.battery.summary.minQgen_W; ...
        v.solar.outflowSummary.min_W;v.solar.spaceRadiationSummary.min_W; ...
        v.solar.structureHeatSummary.min_W;v.solar.totalHeatSummary.min_W;v.solar.storageSummary.min_W];
    Max_W(rows)=[v.pmad.summary.maxQgen_W;v.battery.summary.maxQgen_W; ...
        v.solar.outflowSummary.max_W;v.solar.spaceRadiationSummary.max_W; ...
        v.solar.structureHeatSummary.max_W;v.solar.totalHeatSummary.max_W;v.solar.storageSummary.max_W];
    LoadMin_W(c)=min(v.load_W); LoadMax_W(c)=max(v.load_W);
    MinimumSOC_percent(c)=100*min(v.batterySOC); Unserved_Wh(c)=v.unservedEnergy_Wh;
    Min_C(c)=min(v.solar.temperatureEdges_C); Max_C(c)=max(v.solar.temperatureEdges_C);
    BackMin_C(c)=min(v.solar.backTemperatureEdges_C); BackMax_C(c)=max(v.solar.backTemperatureEdges_C);
    CouplingIterations(c)=v.solarCouplingIterations;
end
tables.heat=table(Case,Quantity,Mean_W,Min_W,Max_W,'VariableNames', ...
    {'Case','Quantity','Mean_W','Min_W','Max_W'});
Case=labels(:);
tables.electrical=table(Case,LoadMin_W,LoadMax_W,MinimumSOC_percent,Unserved_Wh, ...
    'VariableNames',{'Case','LoadMin_W','LoadMax_W','MinimumSOC_percent','Unserved_Wh'});
tables.panelTemperature=table(Case,Min_C,Max_C,BackMin_C,BackMax_C,CouplingIterations,'VariableNames', ...
    {'Case','Min_C','Max_C','BackMin_C','BackMax_C','CouplingIterations'});
Subsystem={'Computer';'ADCS';'Camera';'Active comms'};
Nominal_W=[p.computerNominal_W;p.adcsNominal_W;p.cameraNominal_W;p.commsNominal_W];
Peak_W=[p.computerPeak_W;p.adcsPeak_W;p.cameraPeak_W;p.commsPeak_W];
tables.loads=table(Subsystem,Nominal_W,Peak_W,'VariableNames',{'Subsystem','Nominal_W','Peak_W'});
Setting={'Orbit period';'Propagated history';'Array front area';'Downlinks per orbit'; ...
    'Downlink duration per window';'Camera duty';'Battery series cells';'Battery parallel cells'; ...
    'Effective pack capacity';'Initial SOC';'Minimum allowed SOC';'Maximum allowed SOC'; ...
    'Solar-to-load efficiency';'Charging efficiency';'Battery-to-load efficiency'; ...
    'Assumed cell OCV';'Provisional cell DC resistance';'Pack charge-current limit';'Input battery temperature';'Rear emissivity (assumed)';'Panel thickness (proxy)';'Thermal nodes'};
Value=[r.orbitPeriod_s/60;sum(r.intervalDuration_s)/r.orbitPeriod_s;sp.area_m2; ...
    p.downlinkWindowsPerOrbit;p.downlinkDuration_s/60;100*p.cameraDutyFraction; ...
    bp.seriesCells;bp.parallelCells;r.effectiveCapacity_Ah;100*p.initialSOC; ...
    100*p.minimumSOC;100*p.maximumSOC;100*mp.mpptEfficiency*mp.solarConverterEfficiency; ...
    100*mp.mpptEfficiency*mp.batteryChargeConverterEfficiency;100*mp.batteryDischargeEfficiency; ...
    r.assumedCellOCV_V;bp.cellDCResistance_Ohm;r.packChargeCurrentLimit_A;Tbattery;sp.backEmissivity;sp.panelThickness_m;sp.thermalNodeCount];
Unit={'min';'orbits';'m^2';'windows';'min';'%';'cells';'cells';'Ah';'%';'%';'%'; ...
    '%';'%';'%';'V/cell';'ohm/cell';'A/pack';'deg C';'1';'m';'nodes'};
tables.settings=table(Setting,Value,Unit,'VariableNames',{'Setting','Value','Unit'});
end
function printThermalReport(r)
fprintf('\nCOMPONENT THERMAL SUMMARY\n\nHeat rates [W]\n');
disp(r.summaryTables.heat);
fprintf('Electrical operation\n'); disp(r.summaryTables.electrical);
fprintf('Panel temperature [deg C; Min/Max=front cells; BackMin/Max=rear; includes interval boundaries]\n');
disp(r.summaryTables.panelTemperature);
fprintf('Active component loads [W]\n'); disp(r.summaryTables.loads);
fprintf('Operating settings and assumptions\n'); disp(r.summaryTables.settings);
fprintf(['Table key:\n' ...
    '  PMAD_Qgen / Battery_Qgen  Internal heat generated in the named component.\n' ...
    '  Solar_Qout               Total net heat leaving the array through enabled paths.\n' ...
    '  Solar_Qspace             Net radiation exchanged with space.\n' ...
    '  Solar_QtoStructure       Structure radiation plus mount conduction; positive heats the spacecraft.\n' ...
    '  Solar_QheatInput         Absorbed solar/Earth - terminal electrical export + panel heater.\n' ...
    '  Solar_Qstorage           Rate of thermal energy storage; positive warms the panel, negative cools it.\n' ...
    '  MinimumSOC_percent       Lowest battery charge relative to the effective capacity.\n' ...
    '  Unserved_Wh              Load energy left unmet over the propagated history.\n' ...
    '  CouplingIterations       Iterations used to converge panel temperature and electrical routing.\n\n' ...
    'Notes:\n' ...
    '  QheatInput = Qout + Qstorage. Means are time-weighted over the propagated history.\n' ...
    '  Negative structure transfer heats the panel. Space rejection is not an internal spacecraft heat load.\n']);
enabled=[r.nominal.solar.structureExchangeEnabled,r.peak.solar.structureExchangeEnabled];
if ~any(enabled)
    fprintf('  Structure exchange is DISABLED pending conductance/view factors and sink temperatures; zero is not a verified physical load.\n');
elseif ~all(enabled)
    fprintf('  Structure exchange is enabled for only one operating case; inspect each solar result.\n');
end
fprintf(['  Contacts are evenly spaced assumptions; actual ground-station visibility is not calculated.\n' ...
    '  Blank loads, unknown standby powers and PMAD housekeeping remain unbudgeted.\n' ...
    '  Nominal/peak loads and PMAD efficiencies are provisional. OCV and DC resistance require measured cell data.\n' ...
    '  Front/rear panel temperatures are solved. Battery temperature is an input; battery and PMAD temperatures are not solved.\n' ...
    '  Per-orbit values remain in thermalResult.nominal/peak.perOrbit; tables are in thermalResult.summaryTables.\n\n']);
if r.params.showPerOrbitTables
    for label={'nominal','peak'}
        fprintf('%s per-orbit summary [heat: W; time: min; SOC: %%]\n',label{1});
        disp(struct2table(r.(label{1}).perOrbit));
    end
end
end

%% References and modeling limits
% Operational powers: user-supplied subsystem screenshot; timing: adjustable user-supplied contacts/duty.
% Cell specifications and resistance caveat: BatteryThermalModel.m [BAT], [NASA], [AC/DC].
% PMAD delivered power and losses: PMADThermalModel.m [TI]; solar boundaries: SolarArrayThermalModel.m.
% Discharge-positive circuit: Vterminal=Voc-I*R, Pterminal=I*Voc-I^2*R, dSOC/dt=-I/(3600*capacity_Ah).
% Constant OCV and nominal capacity are placeholders; current/voltage/SOC limits do not constitute a validated charger/BMS design.
% Heat-generation histories are per spacecraft. A constellation-wide sum needs each spacecraft's propagated illumination/modes.
