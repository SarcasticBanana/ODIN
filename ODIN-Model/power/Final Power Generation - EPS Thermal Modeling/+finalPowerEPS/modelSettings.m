function p=modelSettings(group)
% MODELSETTINGS Single editable source of model defaults and assumptions.
% Groups: run, orbit, geometry, array, solar, iv, battery, pmad, component,
% powerBalance, environment, diagnostics. Existing per-call overrides still win.
% Hardware properties, battery assumptions and the 127 C initial state are retained.
% Earth environment constants are published annual-average boundary assumptions.
switch char(group)
case 'array'
geometryDefaults = finalPowerEPS.arrayGeometry(); % shared cell count and provisional CAD coverage

%% Adjustable settings: array studies
% Physical front-face areas and fixed pointing angles; colon notation is start:step:stop.
% Orbit altitude, duration and illumination continue to come from cfg and out.
p.area_m2 = [0.05 0.09 0.10 0.15 0.20];
p.selectedArrayArea_m2 = geometryDefaults.arrayArea_m2;   % project total illuminated face area, not a vendor dimension
p.areaStudyIncidence_deg = 0;       % baseline area study: normal Sun pointing
p.incidenceAngles_deg = 0:1:90;     % each angle is a separate thermal/power case
p.incidenceSweepArea_m2 = p.selectedArrayArea_m2;
p.incidenceSummaryAngles_deg = 0:15:90;
p.runIncidenceSweep = true;        % powerBalance disables the extra study during warm-up

%% Selected product specifications
% [NA] identifies GaInP/GaInAs/Ge cells, coverglass, low-outgassing adhesive and copper planes.
% The wing area, layer thicknesses and assembly thermal properties are not published.
p.productName = 'NanoAvionics 12U double-deployable solar array';
p.productUrl = 'https://nanoavionics.com/cubesat-components/cubesat-gaas-solar-panel/';
p.productRatedWingPower_W = 35;     % per wing; do not treat as the complete-array rating
p.productCellPowerDensity_W_m2 = 368.5; % 36.85 mW/cm^2; area/rating temperature unspecified
p.productPCBThicknessRange_m = [0.0010 0.0016];
p.productOperatingRange_C = [-40 85];

%% Electrical assumptions
% [NA] publishes up to 29.5% cell efficiency; coverage and delivery losses remain estimates.
% The 35 W/wing rating cannot calibrate a physical-area model without the wing's face area.
% Loads and battery behavior are handled separately in powerBalance.m.
p.cellEfficiency = 0.295;           % supplier ceiling, not a guaranteed minimum
p.cellCoverage = geometryDefaults.cellCoverage;              % preserved old CAD coverage; active cell area is unverified
p.referenceTemperature_C = 28;     % retained Spectrolab proxy; [NA] does not give a reference T
p.powerTempCoeff_perC = -0.0021;    % retained proxy, not the selected cell's measured coefficient
p.powerTempCoeffRange_C = [15 75];  % published coefficient range, not operating limits
p.selfShadowRetained = 1.00;
p.mismatchRetained = 0.98;
p.opticalRetained = 0.97;
p.wiringDiodeRetained = 0.97;
p.mpptEfficiency = 0.98;
p.converterEfficiency = 0.95;
p.agingRetained = 0.90;

%% Thermal assumptions (also inherited by SolarArrayThermalModel)
% Two thermal nodes by default: front cells + rear half-plate; direct sunlight heats the front only.
% Optional structure exchange needs actual view factors, temperatures and mounting conductance.
% Radiation to structure assumes large isothermal black surroundings; no enclosure radiosity is solved.
p.initialTemperature_C = 127; % retained hot-start condition, not an average or limit
p.initialBackTemperature_C = 127; % separate rear initial state; warm-up carries both faces forward
p.solarAbsorptivity = 0.91;   % Spectrolab cell value used as a whole-front proxy
p.frontEmissivity = 0.85;     % Spectrolab cell value; exposed substrate finish is unknown
p.spaceTemperature_K = 3;    % ideal cold radiation sink, not a panel property
p.backEmissivity = 0.85;        % exposed rear-face PROXY; replace with measured finish emissivity; 0 disables rear radiation
p.frontStructureViewFactor = 0; % rest of each face views space; 0 omits structure radiation
p.backStructureViewFactor = 0;
p.structureTemperature_C = 20; % placeholder, unused when structure view factors are zero
p.mountConductancePerArea_W_m2K = 0; % omitted until known; enter G_mount/A_front for selected array
p.mountTemperature_C = 20;   % placeholder, unused when mounting conductance is zero
p.absorbedEarthFlux_W_m2 = 0; % optional extra absorbed front-face heat; orbital Earth heating is computed separately

%% Plate property assumptions
% The 1.3 mm carrier is the midpoint of [NA]'s PCB range, not a measured assembly thickness.
% FR-4 is still an unconfirmed laminate proxy; rho, cp and k retain sources [1-4].
% Added layers are excluded until their mass/heat capacity is known; enter them below when available.
p.thermalNodeCount = 2; % 1: uniform plate; 2: front/back half-thickness control volumes
p.panelThickness_m = mean(p.productPCBThicknessRange_m);
p.panelDensity_kg_m3 = 1900;
p.panelSpecificHeat_J_kgK = 1100;
p.panelConductivity_W_mK = 0.25; % through-thickness proxy; node conduction and Biot screening
p.additionalArealMass_kg_m2 = 0; % unknown added cell/glass/copper/adhesive mass, excluded by default
p.additionalArealHeatCapacity_J_m2K = 0; % sum(rho_i*d_i*cp_i) for added layers, once known

%% Solver and diagnostic settings
p.thermalMaxStep_s = 10;
p.thermalRelTol = 1e-7;
p.thermalAbsTol_K = 1e-6;
p.irradianceInterpolation = 'linear';
p.lumpedBiotWarningLimit = 0.1;
p.energyBalanceTolerance_W_m2 = 1e-9;
p.warnTemperatureExtrapolation = true;
p.warnLumpedModel = true;
p.warnProductTemperature = true;

%% Graph and display settings
% makePlots draws figures inside this function; main_demo draws them when makeMainArrayPlots=true.
% Reverse the colormap so normal incidence is red and edge-on incidence is blue.
p.makePlots = false;
p.makeMainArrayPlots = true;
p.showAreaSummary = true;
p.showIncidenceSummary = true;
p.figureColor = [1 1 1];
p.arrayPlotLineWidth = 1.3;
p.mainArrayPlotLineWidth = 1.5;
p.legendLocation = 'southoutside';
p.areaLegendFormat = '%.3g m^2';
p.incidenceFigurePosition = [100 100 1100 600];
p.incidenceMarkerSize = 8;
p.incidenceMarkerAlpha = 0.65;
p.incidenceColormap = 'turbo';
p.incidenceColorLevels = 256;
p.reverseIncidenceColors = true;
p.incidenceColorLimits_deg = [0 90];
p.incidenceColorbarTicks_deg = 0:15:90;

p.environment = finalPowerEPS.modelSettings('environment');
p.frontEarthViewFactor=0; p.backEarthViewFactor=0;
p.frontAbsorbedEarthFlux_W_m2=0; p.backAbsorbedEarthFlux_W_m2=0;

case 'solar'
base=finalPowerEPS.modelSettings('array');
geometryDefaults = finalPowerEPS.arrayGeometry(); % shared cell count and provisional CAD coverage

%% Adjustable array/electrical settings (inherited from arrayPower when supplied)
p.area_m2=geometryDefaults.arrayArea_m2;
p.cellEfficiency=base.cellEfficiency;
p.cellCoverage=base.cellCoverage;
p.referenceTemperature_C=base.referenceTemperature_C;
p.powerTempCoeff_perC=base.powerTempCoeff_perC;
p.powerTempCoeffRange_C=base.powerTempCoeffRange_C;
p.selfShadowRetained=base.selfShadowRetained;
p.mismatchRetained=base.mismatchRetained;
p.opticalRetained=base.opticalRetained;
p.electricalDegradationRetained = 1; % componentThermalOrbit can apply the planning allowance
p.defaultOperatingPointFraction = 1; % used only when actual export is not supplied

%% Adjustable thermal boundaries (defaults match arrayPower)
% Positive heat transfer is NET OUTWARD. Radiation to structure assumes a large black receiver.
% No enclosure radiosity or spacecraft temperature is solved. Boundary values need actual geometry/data.
p.solarAbsorptivity=base.solarAbsorptivity;
p.frontEmissivity=base.frontEmissivity;
p.spaceTemperature_K=base.spaceTemperature_K;
p.backEmissivity=base.backEmissivity;
p.frontStructureViewFactor=base.frontStructureViewFactor;
p.backStructureViewFactor=base.backStructureViewFactor;
p.structureTemperature_C=base.structureTemperature_C;
p.mountConductancePerArea_W_m2K=base.mountConductancePerArea_W_m2K;
p.mountTemperature_C=base.mountTemperature_C;
p.absorbedEarthFlux_W_m2=base.absorbedEarthFlux_W_m2;

%% Adjustable thermal mass/solver (same plate proxies as arrayPower)
p.thermalNodeCount=base.thermalNodeCount;
p.panelThickness_m=base.panelThickness_m;
p.panelDensity_kg_m3=base.panelDensity_kg_m3;
p.panelSpecificHeat_J_kgK=base.panelSpecificHeat_J_kgK;
p.panelConductivity_W_mK=base.panelConductivity_W_mK;
p.lumpedBiotWarningLimit=base.lumpedBiotWarningLimit;
p.additionalArealHeatCapacity_J_m2K=base.additionalArealHeatCapacity_J_m2K;
p.initialTemperature_C=base.initialTemperature_C;
p.initialBackTemperature_C=base.initialBackTemperature_C;
p.thermalMaxStep_s=base.thermalMaxStep_s;
p.energyBalanceTolerance_W_m2=base.energyBalanceTolerance_W_m2;
p.balanceTolerance_W = 1e-7;

%% Adjustable reference case and display
p.defaultNormalIrradiance_W_m2 = 1361;
p.defaultIncidenceAngle_deg = 0;
p.defaultIllumination = 1;
p.defaultCellTemperature_C = 28; % no-history reference point, not predicted temperature
p.defaultHeaterPower_W = 0;
p.productOperatingRange_C=base.productOperatingRange_C;
p.showSummary = true;
p.warnAssumptions = true;

p.buildSummary=true; % omit display-table construction during coupling iterations
p.frontEarthViewFactor=0; p.backEarthViewFactor=0;
p.frontAbsorbedEarthFlux_W_m2=0; p.backAbsorbedEarthFlux_W_m2=0;

case 'battery'
%% Adjustable settings: pack layout and representative load
% The 2S2P arrangement is an editable four-cell design assumption, not a cell specification.
p.seriesCells = 2;
p.parallelCells = 2;
p.defaultBusLoad_W = 23.65;
p.defaultBatteryToBusEfficiency = 0.95; % assumed; used only to infer default current
p.defaultHeaterPower_W = 0;             % unknown; excluded until specified
p.defaultElectronicsPower_W = 0;        % pack-local BMS/balancing heat, not PMAD losses

%% Cell specifications: supplied BatterySpace sheet [BAT]
% Values below are per cell; section 3.2 says 3.7 V although the cover says 3.6 V.
p.cellNominalCapacity_Ah = 2.600;      % [BAT] section 3.1
p.cellMinimumCapacity_Ah = 2.550;      % [BAT] section 3.1
p.cellNominalVoltage_V = 3.7;         % [BAT] section 3.2; cover discrepancy noted above
p.cellStandardChargeCurrent_A = 1.150;% [BAT] sections 3.3 and 6.1, CC/CV
p.cellChargeVoltage_V = 4.2;          % [BAT] section 3.3
p.cellChargeTerminationCurrent_A = 0.050; % [BAT] section 3.3, per cell
p.cellMaxChargeCurrent_A = 2.600;     % [BAT] section 3.4
p.cellStandardDischargeCurrent_A = 0.460; % [BAT] section 3.5
p.cellDischargeCutoffVoltage_V = 3.0; % [BAT] section 3.5
p.cellMaxDischargeCurrent_A = 5.200;  % [BAT] section 3.6, continuous
p.cellMass_kg = 0.0465;               % [BAT] section 3.7, approximate cell-only mass
p.cellChargeTemperatureRange_C = [0 45];     % [BAT] section 3.8
p.cellDischargeTemperatureRange_C = [-20 60];% [BAT] section 3.8
p.cellDiameter_m = 0.0183;            % [BAT] section 4, nominal; maximum 0.0184 m
p.cellHeight_m = 0.0649;              % [BAT] section 4, nominal; maximum 0.0651 m
p.sheetACImpedanceMax_Ohm = 0.090;    % [BAT] section 6.5, <=90 milliohm at 1 kHz, with PTC
p.defaultCellTemperature_C = 23;      % [BAT] section 6.1 test temperature; not a prediction

%% Adjustable thermal assumptions
% DC resistance is NOT supplied by [BAT]; 0.09 ohm is an unvalidated placeholder.
% Its numerical match to the AC limit does not make it a DC value or an upper bound.
p.cellDCResistance_Ohm = 0.090;
p.cellDCResistanceIsMeasured = false;
p.entropyCoefficient_V_K = 0; % d(OCV)/dT unknown; reversible heat omitted by default
p.packInterconnectResistance_Ohm = 0; % unknown pack-local busbar/protection resistance, excluded
p.showSummary = true;
p.warnAssumptions = true;
p.limitTolerance = 1e-10; % numerical comparison tolerance in the stated variable units


case 'pmad'
%% Adjustable settings: provisional path efficiencies
% No PMAD hardware has been selected; these are project assumptions, not vendor ratings.
% Each efficiency covers its full named path; do not add those same losses again.
p.mpptEfficiency = 0.98;
p.solarConverterEfficiency = 0.95;
p.batteryDischargeEfficiency = 0.95;
p.batteryChargeConverterEfficiency = 0.95;
p.extraSwitchResistance_Ohm = 0; % excluded until known; only switches OUTSIDE the path efficiencies

%% Adjustable settings: representative operation
% The default is solar-only load supply with no battery charging or shunt dissipation.
p.defaultSolarToLoad_W = 23.65;
p.defaultBatteryToLoad_W = 0;
p.defaultBatteryCharge_W = 0;
p.defaultControlPower_W = 0; % unknown PMAD housekeeping heat, excluded until specified
p.defaultShuntPower_W = 0;   % unused PV power is NOT automatically dissipated inside PMAD
p.defaultExtraSwitchCurrent_A = 0;
p.showSummary = true;
p.warnAssumptions = true;
p.balanceTolerance_W = 1e-9;


case 'component'
%% Adjustable downlink schedule and camera duty cycle
% These settings apply to BOTH nominal and peak thermal histories.
p.downlinkWindowsPerOrbit = 3; % positive integer; windows are evenly spaced within each orbit
p.downlinkDuration_s = 7.55 * 60; % seconds PER window (5*60 = 5 minutes); 0 disables downlinks
p.cameraDutyFraction = 1.0;    % 0..1 of each orbit: 0=off, 0.25=25%, 1=always on
% Timing options: phases are orbit fractions measured from the propagation start.
p.firstDownlinkPhase = 0; % first window start; 0.5 starts halfway through the first orbit
p.cameraStartPhase = 0;   % camera ON block starts here and lasts cameraDutyFraction*orbitPeriod
p.cameraSunlitOnly = false; % true further restricts camera ON time to illuminated intervals
% Contacts are an assumed schedule, not computed ground-station visibility.
% Camera duty is one continuous ON block per orbit, wrapping at the orbit boundary if needed.
% Computer and ADCS stay on; downlink windows must not overlap.

%% Adjustable operating modes: screenshot supplied by the team
% Lower/upper computer and ADCS values are provisional nominal/peak assignments, not measured mode averages.
% The screenshot gives only one comms/camera value; both cases use it while active. Blank loads remain unbudgeted.
p.computerNominal_W = 3.3;
p.computerPeak_W = 5.0;
p.adcsNominal_W = 1.50;
p.adcsPeak_W = 8.25;
p.cameraNominal_W = 0.6;
p.cameraPeak_W = 0.6;
p.commsNominal_W = 2.0;
p.commsPeak_W = 2.0;
p.commsStandby_W = 0; % unknown off-window consumption, excluded until supplied
p.otherNominal_W = 0; % front end, structures and unspecified loads are NOT yet budgeted
p.otherPeak_W = 0;

%% Adjustable dispatch assumptions
% Solar supplies loads first; a finite battery covers deficits and charges from remaining sunlight.
% SOC is coulomb-counted; constant OCV plus resistance is a preliminary circuit, not a measured OCV-SOC/CC-CV model.
p.initialSOC = 1.0;
p.minimumSOC = 0.40;
p.maximumSOC = 1.0;
p.capacityRetained = 1.0; % BOL capacity; no unverified lifetime derating is imposed
p.cellOpenCircuitVoltage_V = NaN; % NaN uses cell nominal V as an ASSUMED OCV, not a measured OCV
p.batteryCellTemperature_C = NaN; % NaN uses battery model's reference temperature, not predicted temperature
p.packChargeCurrentLimit_A = NaN; % NaN uses Nparallel*datasheet STANDARD charge current
p.externalHarnessRetained = NaN;  % NaN inherits arrayPower wiringDiodeRetained; heat is outside the panel/PMAD
p.usePlanningAgingAllowance = true;% inherit arrayPower agingRetained as an electrical extraction allowance
p.dispatchMaxStep_s = 10;
p.balanceTolerance_W = 1e-7;
p.solarCouplingMaxIterations = 25; % converge electrical routing and panel temperature on the same orbit
p.solarCouplingTolerance_C = 1e-6;
p.solarCouplingTolerance_W = 1e-7;

%% Component settings and display
% Physical component parameters stay at the top of their own function files; optional struct overrides win here.
% The old worst-case powerBalance figure remains a separate study; these cases use out.t/out.illum and out.T.
p.solarOverrides = struct();
p.pmadOverrides = struct();
p.batteryOverrides = struct();
p.showSummary = true;
p.showPerOrbitTables = false;

case 'powerBalance'
%% Adjustable settings: orbit and load
% [BT] gives 35.8 min shadow at 550 km; all penumbra is conservatively treated as dark.
% Combining that shadow with aphelion sunlight is a design envelope, not a propagated date.
p.orbitAltitude_km = 550;
p.eclipseDuration_min = 35.8;
p.sunDistance_AU = 1.0167;          % approximate aphelion distance [ENV]
p.continuousLoad_W = 23.65;        % user-specified: everything on throughout the orbit
p.arrayArea_m2 = [];              % [] inherits arrayPower's selectedArrayArea_m2
p.arrayOverrides = struct();      % optional arrayPower settings; main passes its current p
p.outputStep_s = 5;

%% Adjustable settings: battery
% [BT] Candidate A is an architecture estimate, not a final battery part selection.
% Its 50 W rating is a peak-power screening proxy, not a verified continuous rating.
p.batteryCapacityBOL_Wh = 45;
p.batteryCapacityRetained = 0.80;  % faded-capacity scenario from [BT]
p.batteryInitialSOC = 1.00;        % fully charged at eclipse entry
p.batteryMinimumSOC = 0.40;        % reserve floor relative to faded capacity
p.batteryMaximumSOC = 1.00;
p.batteryMaxDischarge_W = 50;      % pack-side maximum; actual power follows demand
p.batteryToBusEfficiency = 0.95;  % [BT] assumption
p.batteryChargeEfficiency = 1.00; % ideal placeholder: charging hardware is not selected
p.batteryMaxCharge_W = Inf;       % pack-side; unknown limit, ideal unbounded charging
p.allowBatteryInSunlight = false; % requested solar-in-sun / battery-in-shadow split

%% Adjustable settings: thermal settling and display
% Repeat the same thermal orbit until its start/end temperatures agree.
% Only the converged orbit is graphed; the original model's hot-start plots stay unchanged.
p.thermalWarmupMaxOrbits = 20;
p.thermalPeriodicTolerance_C = 0.01;
p.energyTolerance_Wh = 1e-6;
p.showSummary = true;
p.makePlot = true;
p.warnModelLimits = true;
p.figureColor = [1 1 1];
p.figurePosition = [100 100 1100 600];
p.batteryTradeFigureSize_px = [1100 700]; % preferred OUTER size; reduced to fit monitor
p.batteryOperationFigureSize_px = [1100 760];
p.lineWidth = 2.2;
p.eclipseShadeColor = [0.45 0.50 0.60];
p.eclipseShadeAlpha = 0.45;
p.balanceColors = [1.00 0.45 0.00; 0.00 0.40 1.00; 0.80 0.00 0.75; ...
    0.00 0.65 0.30; 0.95 0.05 0.15];
p.balanceLineStyles = {'-', '-', '--', ':', '-.'};
p.balanceLabels = {'Solar only','Battery only','Full-discharge stress test', ...
    'Eclipse support','Extra-load stress test'};
p.supplementalMargin_percent = 10; % target margin; battery limits can prevent reaching it
p.balanceAxisPadding = 1.10;
p.balanceTickStep_percent = 50;

%% Adjustable settings: battery discharge trade
% Each row sets a discharge CEILING; actual discharge follows deficit, solar surplus recharges.
% Highlight the lowest limit meeting interval headroom and restoring initial SOC.
p.runBatteryTrade = true;
p.makeBatteryTradeSurfaceViews = false; % optional old four-view 3D surface; maps/curves are the primary display
p.batteryTradeCurveCommands_percent = [0 25 50 75 100];
p.makeBatteryTradePlots = true; % also requires makePlot = true
p.batteryTradeDischarge_percent = 0:1:100;
p.batteryTradeLineWidth = 1.0;
p.batteryTradeHighlightWidth = 3.0;
p.batteryTradeColorLevels = 256;
p.batteryTradeColorAnchors = [0.00 0.25 1.00; 0.00 0.80 1.00; ...
    1.00 0.85 0.00; 0.95 0.05 0.00]; % blue at zero, red at maximum
p.batteryTradeView_deg = [45 25]; % perspective view; XY, XZ and YZ copies are also opened
p.batteryTradeSurfaceAlpha = 0.75;
p.batteryTradeMeshRateStep_percent = 10;
p.batteryTradeMeshTimeStep_min = 5;
p.batteryTradeMeshLineWidth = 0.45;
p.batteryTradeMeshColor = [0.25 0.25 0.25];
p.batteryTradeTolerance_percent = 1e-8;


case 'iv'
geometry=finalPowerEPS.arrayGeometry();
%% Wiring assumptions and geometry (Ns*Np must equal totalCells)
p.mode='power_consistent_proxy';
p.totalCells=geometry.totalCells;
p.seriesCells=16; % ASSUMED scenario, not inferred from CAD or confirmed wiring
p.parallelStrings=[]; % [] = totalCells/seriesCells; explicitly enter a count only if known
p.arrayArea_m2=geometry.arrayArea_m2; p.cellCoverage=geometry.cellCoverage;
p.testVoltage_V=[]; % optional operating-point test; no fixed voltage by default
p.curvePoints=301; p.makePlots=true; p.showSummary=true; p.snapshotCount=4;
% Display only: full outer window is capped to the current monitor by the helper.
p.figureSize_px=[1120 740]; p.plotLineWidth=1.9; p.plotFontSize=10;
p.plotPaddingFraction=0.12; % space above maxima; includes the whole current history

%% Spectrolab XTJ Prime BOL reference cell, manufacturer sheet (2018)
% https://www.spectrolab.com/photovoltaics/XTJ-Prime_Data_Sheet.pdf
p.referenceIrradiance_W_m2=1353; p.referenceTemperature_C=28;
p.cellVoc_V=2.720; p.cellVmp_V=2.406;
p.cellJsc_A_cm2=.018; p.cellJmp_A_cm2=.0175;
p.cellVocTemperature_V_K=-.0056; p.cellJscTemperature_A_cm2K=10e-6;
p.cellVmpTemperature_V_K=-.0063; p.cellJmpTemperature_A_cm2K=5e-6;
p.coefficientRange_C=[15 125];
% MPP temperature coefficients are retained comparison values, not exact extra fit constraints.

case 'orbit'
p.mu = 398600.4418; % p.mu: Earth gravitational parameter [km^3/s^2]
p.Re = 6378.137; % p.Re: Earth equatorial radius [km]
p.J2 = 1.082626683e-3; % p.J2: Earth second zonal harmonic coefficient
p.Rsun = 695700; % p.Rsun: Sun mean radius [km]
p.AU = 149597870.7; % p.AU: astronomical unit [km]
p.S0 = 1361; % p.S0: solar constant, irradiance at 1 AU [W/m^2]

p.epoch = [2026 9 15 0 0 0]; % p.epoch: UTC calendar instant at sim start [yr mon day hr min s]

p.h = 800; % p.h: orbit altitude above p.Re [km]
p.inc = 98.6; % p.inc: orbit inclination [deg]
p.raan0 = 0; % p.raan0: initial right ascension of ascending node [deg]
p.argp0 = 0; % p.argp0: initial argument of perigee [deg]
p.nu0 = 0; % p.nu0: initial true anomaly [deg]
p.ecc = 0; % p.ecc: orbit eccentricity

p.use.J2 = true; % p.use.J2: enable the J2 acceleration term
p.use.eclipse = true; % p.use.eclipse: enable the Sun direction and Earth shadow history
p.use.drag = false; % p.use.drag: enable atmospheric drag [not implemented]
p.use.srp = false; % p.use.srp: enable solar radiation pressure [not implemented]

p.eclipse.model = 'conical'; % p.eclipse.model: shadow geometry, 'conical' or 'cylindrical'

p.ode.relTol = 1e-11; % p.ode.relTol: ode45 relative error tolerance
p.ode.absTol = 1e-11; % p.ode.absTol: ode45 absolute error tolerance
p.ode.dtOut = 60; % p.ode.dtOut: spacing of the uniform output time grid [s]

p.eclipse.searchStep_s=10; % bracket shadow events independently of output plotting cadence
p.eclipse.transitionStep_s=1; % resolve partial illumination between refined contacts

case 'geometry'
p.totalCells=64;
p.legacyArrayArea_m2=0.1999;
p.legacyCellArea_m2=0.00310975;
p.defaultArrayArea_m2=0.175;

case 'environment'
% Uniform diffuse spherical Earth. Annual averages, not worst-case design bounds.
% NASA SP-8067: https://ntrs.nasa.gov/citations/19710023628
p.enabled=true;
p.earthInfrared_W_m2=237;
p.earthAlbedo=0.30;
p.radialOrder=16; p.azimuthCount=64;
% Without out.panelNormal_ECI, use the original Sun-facing array assumption.
% Fixed-angle studies tilt toward the orbit normal (an explicit azimuth convention).
p.normalMode='sun_tracking';
case 'diagnostics'
p.makeAdditionalPlots=true;
p.figureSize_px=[1120 760];
case 'run'
p.altitude_km=600;
p.numberOfOrbits=10;
otherwise
error('modelSettings:group','Unknown settings group: %s',group);
end
end
