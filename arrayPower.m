function result = arrayPower(out, cfg, area_m2, overrides)
% ARRAYPOWER Area sizing and incidence sweeps with a one-sided transient plate.
% Areas are physical front-face areas [m^2]; the baseline defaults to normal Sun incidence.
% Requires solarFlux.m and the existing orbital-model output.

%% Usage
% Edit the settings below, then run main_demo.m or result = arrayPower(out, cfg).
% Optional third argument replaces the default area list; fourth argument overrides p fields.
% result = arrayPower(out, cfg, [0.05 0.10 0.15 0.20]);
% changes = struct('initialTemperature_C', 127, 'makePlots', true);
% result = arrayPower(out, cfg, [0.05 0.10], changes);
% disp(result.summary)
% disp(result.incidenceSweep.summary)

%% Adjustable settings: array studies
% Physical front-face areas and fixed pointing angles; colon notation is start:step:stop.
% Orbit altitude, duration and illumination continue to come from cfg and out.
p.area_m2 = [0.05 0.09 0.10 0.15 0.20];
p.areaStudyIncidence_deg = 0;       % baseline area study: normal Sun pointing
p.incidenceAngles_deg = 0:1:90;     % each angle is a separate thermal/power case
p.incidenceSweepArea_m2 = 0.199;
p.incidenceSummaryAngles_deg = 0:15:90;

%% Electrical assumptions
% The default area study retains the electrical factors and normal Sun pointing.
% Planning output retains the same 10% allowance, with no assigned mission life.
% Loads, battery behavior and EPS capacity limits are still excluded.
p.cellEfficiency = 0.30;
p.cellCoverage = 0.80;
p.referenceTemperature_C = 28;
p.powerTempCoeff_perC = -0.0021;
p.powerTempCoeffRange_C = [15 75];  % published coefficient range, not operating limits
p.selfShadowRetained = 1.00;
p.mismatchRetained = 0.98;
p.opticalRetained = 0.97;
p.wiringDiodeRetained = 0.97;
p.mpptEfficiency = 0.98;
p.converterEfficiency = 0.95;
p.agingRetained = 0.90;

%% Thermal assumptions
% One uniform plate temperature; only the illuminated front face can radiate.
% That face continues radiating in eclipse, with unit view factor to space.
% The back/edges are adiabatic; Earth heating, mounting conduction and convection are omitted.
p.initialTemperature_C = 127; % retained hot-start condition, not an average or limit
p.solarAbsorptivity = 0.91;   % Spectrolab cell value used as a whole-front proxy
p.frontEmissivity = 0.85;     % Spectrolab cell value; exposed substrate finish is unknown
p.spaceTemperature_K = 3;    % ideal cold radiation sink, not a panel property

%% Plate property assumptions
% Provisional FR-4 carrier plate; cells, coverglass, adhesives and stiffeners add thermal mass.
% Density, cp and k are rounded means of four published FR-4 model inputs [1-4].
% The 1.6 mm carrier thickness [5-6] excludes those added layers; replace with assembly data later.
p.panelThickness_m = 0.0016;
p.panelDensity_kg_m3 = 1900;
p.panelSpecificHeat_J_kgK = 1100;
p.panelConductivity_W_mK = 0.25; % through-thickness proxy; lumped-model check only

%% Solver and diagnostic settings
p.thermalMaxStep_s = 10;
p.thermalRelTol = 1e-7;
p.thermalAbsTol_K = 1e-6;
p.irradianceInterpolation = 'linear';
p.lumpedBiotWarningLimit = 0.1;
p.energyBalanceTolerance_W_m2 = 1e-9;
p.warnTemperatureExtrapolation = true;
p.warnLumpedModel = true;

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
p.legendLocation = 'best';
p.areaLegendFormat = '%.3g m^2';
p.incidenceFigurePosition = [100 100 1100 600];
p.incidenceMarkerSize = 8;
p.incidenceMarkerAlpha = 0.65;
p.incidenceColormap = 'turbo';
p.incidenceColorLevels = 256;
p.reverseIncidenceColors = true;
p.incidenceColorLimits_deg = [0 90];
p.incidenceColorbarTicks_deg = 0:15:90;

%% Physical constant
% Unit conversions and mathematical constants in the equations are fixed definitions.
sigma = 5.670374419e-8;        % Stefan-Boltzmann constant [W/(m^2 K^4)]

%% Parameter checks
% Apply optional overrides, then check scalar, vector and text settings separately.
if nargin < 3
    area_m2 = [];
end
if nargin < 4
    overrides = struct();
    if isstruct(area_m2)
        overrides = area_m2;
        area_m2 = [];
    end
end
validateattributes(overrides, {'struct'}, {'scalar'});
names = fieldnames(overrides);
for k = 1:numel(names)
    if ~isfield(p, names{k})
        error('arrayPower:parameter', 'Unknown parameter: %s', names{k});
    end
    p.(names{k}) = overrides.(names{k});
end
if ~isempty(area_m2)
    p.area_m2 = area_m2;
end
vectorFields = {'area_m2','incidenceAngles_deg','incidenceSummaryAngles_deg', ...
    'powerTempCoeffRange_C','figureColor','incidenceFigurePosition', ...
    'incidenceColorLimits_deg','incidenceColorbarTicks_deg'};
textFields = {'irradianceInterpolation','legendLocation','areaLegendFormat', ...
    'incidenceColormap'};
names = fieldnames(p);
for k = 1:numel(names)
    name = names{k};
    if ismember(name, textFields)
        if ~(ischar(p.(name)) && isrow(p.(name))) && ...
                ~(isstring(p.(name)) && isscalar(p.(name)))
            error('arrayPower:parameter', '%s must be a text scalar.', name);
        end
        p.(name) = char(p.(name));
        validateattributes(p.(name), {'char'}, {'nonempty','row'});
    elseif ismember(name, vectorFields)
        validateattributes(p.(name), {'numeric'}, {'real','finite','vector','nonempty'});
    else
        validateattributes(p.(name), {'numeric','logical'}, {'real','finite','scalar'});
    end
end
fractionFields = {'cellEfficiency','cellCoverage','selfShadowRetained', ...
    'mismatchRetained','opticalRetained','wiringDiodeRetained', ...
    'mpptEfficiency','converterEfficiency','agingRetained', ...
    'solarAbsorptivity','frontEmissivity','incidenceMarkerAlpha'};
for k = 1:numel(fractionFields)
    validateattributes(double(p.(fractionFields{k})), {'numeric'}, ...
        {'>=',0,'<=',1});
end
positiveFields = {'panelThickness_m','panelDensity_kg_m3', ...
    'panelSpecificHeat_J_kgK','panelConductivity_W_mK', ...
    'thermalMaxStep_s','thermalRelTol','thermalAbsTol_K','lumpedBiotWarningLimit', ...
    'arrayPlotLineWidth','mainArrayPlotLineWidth','incidenceMarkerSize'};
for k = 1:numel(positiveFields)
    validateattributes(p.(positiveFields{k}), {'numeric'}, {'positive'});
end
validateattributes(p.initialTemperature_C, {'numeric'}, {'>',-273.15});
validateattributes(p.referenceTemperature_C, {'numeric'}, {'>',-273.15});
validateattributes(p.spaceTemperature_K, {'numeric'}, {'nonnegative'});
% Numeric bounds require numeric values, including when a flag is logical.
flagFields = {'makePlots','makeMainArrayPlots','showAreaSummary','showIncidenceSummary', ...
    'reverseIncidenceColors','warnTemperatureExtrapolation','warnLumpedModel'};
for k = 1:numel(flagFields)
    validateattributes(double(p.(flagFields{k})), {'numeric'}, {'integer','>=',0,'<=',1});
    p.(flagFields{k}) = logical(p.(flagFields{k}));
end
validateattributes(p.area_m2, {'numeric'}, {'nonnegative'});
validateattributes(p.incidenceSweepArea_m2, {'numeric'}, {'nonnegative'});
validateattributes(p.energyBalanceTolerance_W_m2, {'numeric'}, {'nonnegative'});
validateattributes(p.areaStudyIncidence_deg, {'numeric'}, {'>=',0,'<=',90});
angleFields = {'incidenceAngles_deg','incidenceSummaryAngles_deg', ...
    'incidenceColorLimits_deg','incidenceColorbarTicks_deg'};
for k = 1:numel(angleFields)
    validateattributes(p.(angleFields{k}), {'numeric'}, {'>=',0,'<=',90});
end
validateattributes(p.powerTempCoeffRange_C, {'numeric'}, {'numel',2,'>',-273.15});
validateattributes(p.incidenceColorLimits_deg, {'numeric'}, {'numel',2});
if diff(p.powerTempCoeffRange_C) <= 0 || diff(p.incidenceColorLimits_deg) <= 0
    error('arrayPower:parameter', 'Temperature range and color limits must increase.');
end
validateattributes(p.figureColor, {'numeric'}, {'numel',3,'>=',0,'<=',1});
validateattributes(p.incidenceFigurePosition, {'numeric'}, {'numel',4});
validateattributes(p.incidenceFigurePosition(3:4), {'numeric'}, {'positive'});
validateattributes(p.incidenceColorLevels, {'numeric'}, {'integer','>=',2});
p.irradianceInterpolation = validatestring(p.irradianceInterpolation, ...
    {'linear','nearest','previous','next','pchip'});
% Plot properties use row vectors even if an override was supplied as a column.
p.figureColor = p.figureColor(:).';
p.incidenceFigurePosition = p.incidenceFigurePosition(:).';
p.incidenceColorLimits_deg = p.incidenceColorLimits_deg(:).';
p.incidenceColorbarTicks_deg = p.incidenceColorbarTicks_deg(:).';

%% Orbital inputs
% solarFlux already includes Earth eclipse and Sun distance; do not apply eclipse twice.
if ~isstruct(out) || ~isscalar(out) || ~isfield(out,'t') || ~isfield(out,'illum')
    error('arrayPower:input', 'out must contain time and illumination histories.');
end
t = out.t(:).';
illum = out.illum(:).';
validateattributes(t, {'numeric'}, {'real','finite','vector','nonempty'});
validateattributes(illum, {'numeric'}, ...
    {'real','finite','vector','>=',0,'<=',1});
if numel(t) < 2 || any(diff(t) <= 0) || numel(illum) ~= numel(t)
    error('arrayPower:time', 'Need >=2 increasing times and matching illumination.');
end
t = double(t);
illum = double(illum);
S = solarFlux(out, cfg);
validateattributes(S, {'numeric'}, {'real','finite','vector','nonnegative'});
S = double(S(:).');
if numel(S) ~= numel(t)
    error('arrayPower:flux', 'Irradiance and time histories must have equal lengths.');
end
Snormal = S;
S = Snormal * cosd(p.areaStudyIncidence_deg);

%% Transient heat section
% C_A*dT/dt = alpha*S_front - P_extracted/A - epsilon*sigma*(T^4-T_space^4).
% All radiation temperatures are in kelvin; C_A = rho*thickness*cp [J/(m^2 K)].
% Identical plate construction makes area cancel, so all area cases share T(t).
arealMass = p.panelDensity_kg_m3 * p.panelThickness_m;
arealHeatCapacity = arealMass * p.panelSpecificHeat_J_kgK;
elapsed = t - t(1);
opts = odeset('RelTol', p.thermalRelTol, 'AbsTol', p.thermalAbsTol_K, ...
    'MaxStep', min(p.thermalMaxStep_s, min(diff(elapsed))));
rhs = @(time,T) plateHeatRate(time, T, elapsed, S, p, sigma, arealHeatCapacity);
solution = ode45(rhs, [0 elapsed(end)], p.initialTemperature_C + 273.15, opts);
if solution.x(end) < elapsed(end)
    error('arrayPower:thermalSolver', 'Thermal solver did not reach the window end.');
end
temperature_K = deval(solution, elapsed);
if any(~isfinite(temperature_K)) || any(temperature_K <= 0)
    error('arrayPower:temperature', 'Thermal solution must remain finite and above 0 K.');
end
cellTemperature_C = temperature_K - 273.15;

%% Power section
% Electrical output follows the solved temperature without a hot/cold switch.
% Planning remains a multiplier on BOL output, with no separate aged thermal solution.
arrayBOL_per_m2 = arrayPowerDensity(S, cellTemperature_C, p);
busBOL_per_m2 = arrayBOL_per_m2 * p.wiringDiodeRetained * ...
    p.mpptEfficiency * p.converterEfficiency;
isGenerating = S > 0 & p.selfShadowRetained > 0 & ...
    p.cellCoverage > 0 & p.cellEfficiency > 0;
if p.warnTemperatureExtrapolation && any(isGenerating & ...
        (cellTemperature_C < p.powerTempCoeffRange_C(1) | ...
         cellTemperature_C > p.powerTempCoeffRange_C(2)))
    warning('arrayPower:temperature', ...
        ['Generating temperatures extend outside the coefficient range ' ...
         '%g to %g C; the electrical correction is extrapolated.'], ...
         p.powerTempCoeffRange_C(1), p.powerTempCoeffRange_C(2));
end

%% Thermal checks section
% The radiation Biot number screens the equivalent plate, not adhesive/contact layers.
% Lc = volume/radiating area = thickness because only one face radiates.
hRadiation = p.frontEmissivity * sigma * ...
    (temperature_K + p.spaceTemperature_K) .* ...
    (temperature_K.^2 + p.spaceTemperature_K^2);
biot = hRadiation * p.panelThickness_m / p.panelConductivity_W_mK;
if p.warnLumpedModel && any(biot > p.lumpedBiotWarningLimit)
    warning('arrayPower:lumpedModel', ...
        'Radiation Biot number exceeds %g; check a layered conduction model.', ...
        p.lumpedBiotWarningLimit);
end
absorbedSolar = p.solarAbsorptivity * S * p.selfShadowRetained;
extractedElectrical = arrayBOL_per_m2 * p.mpptEfficiency;
netRadiation = p.frontEmissivity * sigma * ...
    (temperature_K.^4 - p.spaceTemperature_K^4);
netHeat = absorbedSolar - extractedElectrical - netRadiation;

%% Results section
% Existing power fields retain their names and area-by-time matrix dimensions.
% Thermal fluxes are per square metre of physical panel face.
A = double(p.area_m2(:));
result.time_s = t;
result.area_m2 = A;
result.irradiance_W_m2 = S;
result.normalIrradiance_W_m2 = Snormal;
result.cellTemperature_C = cellTemperature_C;
result.arrayBOL_W = A * arrayBOL_per_m2;
result.arrayPlanning_W = result.arrayBOL_W * p.agingRetained;
result.busBOL_W = A * busBOL_per_m2;
result.busPlanning_W = result.busBOL_W * p.agingRetained;
result.params = p;
result.windowDuration_s = elapsed(end);
result.equivalentSunFraction = trapz(t, illum) / result.windowDuration_s;
result.thermal.arealMass_kg_m2 = arealMass;
result.thermal.arealHeatCapacity_J_m2K = arealHeatCapacity;
result.thermal.absorbedSolar_W_m2 = absorbedSolar;
result.thermal.extractedElectricalBOL_W_m2 = extractedElectrical;
result.thermal.netRadiation_W_m2 = netRadiation;
result.thermal.netHeat_W_m2 = netHeat;
result.thermal.temperatureRate_K_s = netHeat / arealHeatCapacity;
result.thermal.radiationBiot = biot;
result.thermal.maxRadiationBiot = max(biot);
result.thermal.solarThroughThicknessScale_K = ...
    max(absorbedSolar) * p.panelThickness_m / p.panelConductivity_W_mK;
result.thermal.temperatureMin_C = min(cellTemperature_C);
result.thermal.temperatureMax_C = max(cellTemperature_C);

%% Energy summary section
% Energy and average power cover the supplied time window, not necessarily one orbit.
energyBOL_Wh = trapz(t, result.busBOL_W, 2) / 3600;
energyPlanning_Wh = trapz(t, result.busPlanning_W, 2) / 3600;
averageBOL_W = energyBOL_Wh * 3600 / result.windowDuration_s;
averagePlanning_W = energyPlanning_Wh * 3600 / result.windowDuration_s;
result.summary = table(A, max(result.busBOL_W,[],2), ...
    max(result.busPlanning_W,[],2), averageBOL_W, averagePlanning_W, ...
    energyBOL_Wh, energyPlanning_Wh, 'VariableNames', ...
    {'Area_m2','PeakBOL_W','PeakPlanning_W','WindowAvgBOL_W', ...
     'WindowAvgPlanning_W','WindowEnergyBOL_Wh','WindowEnergyPlanning_Wh'});

%% Incidence sweep section
% Each angle stays fixed relative to the Sun and has its own temperature history.
% Both studies take their pointing and area settings from the parameter section.
result.incidenceSweep = solveIncidenceSweep(t, Snormal, p, sigma, arealHeatCapacity, opts);

%% Graphing section
% Set makePlots=true to display temperature and planning bus power in separate figures.
if p.makePlots
    figure('Color',p.figureColor,'Name','Array temperature');
    plot(t/60, cellTemperature_C, 'LineWidth', p.arrayPlotLineWidth);
    grid on; xlabel('Time [min]'); ylabel('Cell temperature [deg C]');
    title('One-sided transient array temperature');
    figure('Color',p.figureColor,'Name','Array power');
    plot(t/60, result.busPlanning_W.', 'LineWidth', p.arrayPlotLineWidth);
    grid on; xlabel('Time [min]'); ylabel('Planning bus power [W]');
    title('Array power with transient temperature and aging allowance');
    labels = arrayfun(@(a) sprintf(p.areaLegendFormat, a), A, 'UniformOutput', false);
    legend(labels, 'Location',p.legendLocation);
end

%% Incidence graphing section
% Add power and temperature figures; every point is colored by its incidence angle.
if p.makePlots
    plotIncidenceSweep(result.incidenceSweep, p);
end
end

%% Heat-balance function
% Interpolate the supplied irradiance with the selected method; missed eclipses require finer inputs.
% MPPT shortfall stays as panel heat; harness/converter dissipation is assumed outside the plate.
function dTdt = plateHeatRate(time, T_K, t, S, p, sigma, heatCapacity, incidenceFactor)
% Only the seven-input baseline call omits the incidence factor.
if nargin < 8
    incidenceFactor = 1;
end
time = min(max(time, t(1)), t(end));
flux = interp1(t, S, time, p.irradianceInterpolation) .* incidenceFactor;
absorbed = p.solarAbsorptivity * flux * p.selfShadowRetained;
extracted = arrayPowerDensity(flux, T_K - 273.15, p) * p.mpptEfficiency;
radiated = p.frontEmissivity * sigma * (T_K.^4 - p.spaceTemperature_K^4);
dTdt = (absorbed - extracted - radiated) / heatCapacity;
end

%% Electrical power function
% Use the existing temperature coefficient only where electrical generation is possible.
% Absorptivity already represents whole-face solar absorption; do not multiply it by opticalRetained.
function available_W_m2 = arrayPowerDensity(S, temperature_C, p)
active = S > 0 & p.selfShadowRetained > 0 & ...
    p.cellCoverage > 0 & p.cellEfficiency > 0;
temperatureFactor = ones(size(S));
temperatureFactor(active) = 1 + p.powerTempCoeff_perC * ...
    (temperature_C(active) - p.referenceTemperature_C);
cellEfficiency = p.cellEfficiency * temperatureFactor;
if any(cellEfficiency(active) < 0 | cellEfficiency(active) > 1)
    error('arrayPower:efficiency', ...
        'Extrapolated cell efficiency is outside [0,1]; revise the temperature model/data.');
end
available_W_m2 = S .* cellEfficiency * p.cellCoverage * ...
    p.selfShadowRetained * p.mismatchRetained * p.opticalRetained;
absorbed_W_m2 = p.solarAbsorptivity * S * p.selfShadowRetained;
if any(available_W_m2 > absorbed_W_m2 + p.energyBalanceTolerance_W_m2, 'all')
    error('arrayPower:energyBalance', ...
        'Electrical output exceeds absorbed sunlight; check efficiency and absorptivity.');
end
end

%% Incidence thermal function
% theta is measured from the outward front normal to the direction toward the Sun.
% Photon travel is opposite that Sun direction; theta=0 faces the Sun and 90 is edge-on.
% Projected sunlight drives heat and power, while radiating area and thermal mass stay fixed.
function sweep = solveIncidenceSweep(t, S, p, sigma, heatCapacity, opts)
sweep.area_m2 = double(p.incidenceSweepArea_m2);
sweep.angle_deg = double(p.incidenceAngles_deg(:));
sweep.cosineFactor = cosd(sweep.angle_deg);
sweep.time_s = t;
sweep.irradiance_W_m2 = sweep.cosineFactor * S;
elapsed = t - t(1);
sweep.windowDuration_s = elapsed(end);

% Solve independent angle cases together using the same heat balance as the area study.
rhs = @(time,T) plateHeatRate(time, T, elapsed, S, p, sigma, ...
    heatCapacity, sweep.cosineFactor);
initial_K = repmat(p.initialTemperature_C + 273.15, numel(sweep.angle_deg), 1);
solution = ode45(rhs, [0 elapsed(end)], initial_K, opts);
if solution.x(end) < elapsed(end)
    error('arrayPower:incidenceSolver', 'Incidence solver did not reach the window end.');
end
temperature_K = deval(solution, elapsed);
if any(~isfinite(temperature_K) | temperature_K <= 0, 'all')
    error('arrayPower:incidenceTemperature', 'Incidence temperatures must be finite and above 0 K.');
end
sweep.cellTemperature_C = temperature_K - 273.15;

%% Incidence power section
% Every angle uses its own projected irradiance and solved cell temperature.
% Output matrices are angle-by-time; the planning allowance matches the original study.
available = arrayPowerDensity(sweep.irradiance_W_m2, sweep.cellTemperature_C, p);
sweep.arrayBOL_W = sweep.area_m2 * available;
sweep.arrayPlanning_W = sweep.arrayBOL_W * p.agingRetained;
sweep.busBOL_W = sweep.arrayBOL_W * p.wiringDiodeRetained * ...
    p.mpptEfficiency * p.converterEfficiency;
sweep.busPlanning_W = sweep.busBOL_W * p.agingRetained;
active = sweep.irradiance_W_m2 > 0 & p.selfShadowRetained > 0 & ...
    p.cellCoverage > 0 & p.cellEfficiency > 0;
sweep.temperatureCoefficientExtrapolated = any(active & ...
    (sweep.cellTemperature_C < p.powerTempCoeffRange_C(1) | ...
     sweep.cellTemperature_C > p.powerTempCoeffRange_C(2)), 2);
if p.warnTemperatureExtrapolation && any(sweep.temperatureCoefficientExtrapolated)
    warning('arrayPower:incidenceExtrapolation', ...
        ['Some incidence cases generate outside the %g to %g C coefficient range; ' ...
         'their electrical temperature correction is extrapolated.'], ...
         p.powerTempCoeffRange_C(1), p.powerTempCoeffRange_C(2));
end

%% Incidence thermal checks section
% Heat fluxes use physical face area; radiation has no solar-incidence cosine factor.
sweep.thermal.arealHeatCapacity_J_m2K = heatCapacity;
sweep.thermal.absorbedSolar_W_m2 = p.solarAbsorptivity * ...
    sweep.irradiance_W_m2 * p.selfShadowRetained;
sweep.thermal.extractedElectricalBOL_W_m2 = available * p.mpptEfficiency;
sweep.thermal.netRadiation_W_m2 = p.frontEmissivity * sigma * ...
    (temperature_K.^4 - p.spaceTemperature_K^4);
sweep.thermal.netHeat_W_m2 = sweep.thermal.absorbedSolar_W_m2 - ...
    sweep.thermal.extractedElectricalBOL_W_m2 - sweep.thermal.netRadiation_W_m2;
sweep.thermal.temperatureRate_K_s = sweep.thermal.netHeat_W_m2 / heatCapacity;
hRadiation = p.frontEmissivity * sigma * (temperature_K + p.spaceTemperature_K) .* ...
    (temperature_K.^2 + p.spaceTemperature_K^2);
sweep.thermal.maxRadiationBiot = max(hRadiation, [], 2) * ...
    p.panelThickness_m / p.panelConductivity_W_mK;
if p.warnLumpedModel && any(sweep.thermal.maxRadiationBiot > p.lumpedBiotWarningLimit)
    warning('arrayPower:incidenceLumpedModel', ...
        'An incidence case exceeds Bi=%g; check a layered conduction model.', ...
        p.lumpedBiotWarningLimit);
end

%% Incidence summary section
% Window energy includes the shared initial temperature and every supplied orbit sample.
energyBOL_Wh = trapz(t, sweep.busBOL_W, 2) / 3600;
energyPlanning_Wh = trapz(t, sweep.busPlanning_W, 2) / 3600;
averagePlanning_W = energyPlanning_Wh * 3600 / sweep.windowDuration_s;
sweep.summary = table(sweep.angle_deg, max(sweep.busPlanning_W, [], 2), ...
    averagePlanning_W, energyBOL_Wh, energyPlanning_Wh, ...
    min(sweep.cellTemperature_C, [], 2), max(sweep.cellTemperature_C, [], 2), ...
    sweep.temperatureCoefficientExtrapolated, 'VariableNames', ...
    {'Angle_deg','PeakPlanning_W','WindowAvgPlanning_W','WindowEnergyBOL_Wh', ...
     'WindowEnergyPlanning_Wh','MinTemperature_C','MaxTemperature_C', ...
     'TempCoeffExtrapolated'});
[~, bestIndex] = max(energyPlanning_Wh);
sweep.bestEnergyAngle_deg = sweep.angle_deg(bestIndex);
end

%% Incidence plotting function
% Both figures use all angle/time samples and the same angle color scale.
function plotIncidenceSweep(sweep, p)
[timeGrid_min, angleGrid_deg] = meshgrid(sweep.time_s / 60, sweep.angle_deg);
colors = feval(p.incidenceColormap, p.incidenceColorLevels);
if p.reverseIncidenceColors
    colors = flipud(colors);
end
figure('Color',p.figureColor,'Name','Array incidence power', ...
    'Position',p.incidenceFigurePosition);
scatter(timeGrid_min(:), sweep.busPlanning_W(:), p.incidenceMarkerSize, angleGrid_deg(:), ...
    'filled', 'MarkerFaceAlpha', p.incidenceMarkerAlpha);
grid on; xlabel('Time [min]'); ylabel('Planning bus power [W]');
title(sprintf('Coupled incidence sweep: power, A = %g m^2', sweep.area_m2));
colormap(gca, colors); caxis(p.incidenceColorLimits_deg);
cb = colorbar; cb.Label.String = 'Incidence angle [deg]'; cb.Ticks = p.incidenceColorbarTicks_deg;

figure('Color',p.figureColor,'Name','Array incidence temperature', ...
    'Position',p.incidenceFigurePosition);
scatter(timeGrid_min(:), sweep.cellTemperature_C(:), p.incidenceMarkerSize, angleGrid_deg(:), ...
    'filled', 'MarkerFaceAlpha', p.incidenceMarkerAlpha);
grid on; xlabel('Time [min]'); ylabel('Cell temperature [deg C]');
title(sprintf('Coupled incidence sweep: temperature, A = %g m^2', sweep.area_m2));
colormap(gca, colors); caxis(p.incidenceColorLimits_deg);
cb = colorbar; cb.Label.String = 'Incidence angle [deg]'; cb.Ticks = p.incidenceColorbarTicks_deg;
end

%% References
% Sources checked 2026-10-01; [1-6] support the carrier-property assumptions above.
% Each entry identifies the quantities used; project choices are listed separately below.

%% Heat-transfer references
% [HT] Cengel, Y. A., Cimbala, J. M., and Turner, R. H.
% Fundamentals of Thermal-Fluid Sciences, 5th ed., McGraw-Hill Education.
% Uploaded Ch. 18: Eq. 18-1, p. 718 (lumped energy storage);
% Eq. 18-9 and Bi <= 0.1 criterion, pp. 720-721 (lumped-model screening).
% Uploaded Ch. 21: Eq. 21-53 in Table 21-5, p. 897 (radiation to large surroundings).
% Publisher's edition/author listing:
% https://www.mheducation.com.au/ebook-fundamental-of-thermal-fluid-sciences-in-si-units-5th-edition-9789814923194-aus

%% Heat-balance derivation
% The model applies [HT]'s energy balance with solar input, electrical removal and radiation:
% rho*d*cp*dT/dt = alpha*S_front - P_extracted/A - epsilon*sigma*(T^4-T_space^4).
% This nonlinear ODE is solved numerically; the chapter's convection-only exponential is not used.
%
% Factoring T^4-T_space^4 gives h_rad = epsilon*sigma*(T+T_space)*(T^2+T_space^2).
% With one radiating face, Lc = V/A = d; Bi = h_rad*d/k is a screening estimate.
% The reported alpha*S*d/k temperature scale follows 1D Fourier conduction, not a layer solution.

%% Solar-cell reference
% [CELL] Spectrolab, Inc. "32% XTE+ LEO Space Qualified Triple Junction Solar Cell."
% Data sheet, pp. 1-2: T_ref=28 C, BOL efficiency=0.322 at AM0=135.3 mW/cm^2;
% BOL dPmp/dT=-92 microW/(cm^2 C), fitted over 15-75 C; alpha=0.91, epsilon=0.85.
% https://www.spectrolab.com/photovoltaics/XTE%2B%20LEO%20Data%20Sheet.pdf
%
% gamma = (dPmp/dT)/Pmp_ref = -92e-6/(0.322*135.3e-3) = -0.00211171/C.
% The code rounds gamma to -0.0021/C and uses eta(T)=eta_ref*[1+gamma*(T_C-28)].
% Applying cell optics to the whole front and this slope to a 30% cell are project approximations.
%
% The 15-75 C range applies to the published coefficient, not a panel operating limit.
% Results outside that range extrapolate the linear correction; no aged coefficient is used.

%% Thermal-property sources
% FR-4 values below are rho [kg/m^3], cp [J/(kg K)], k [W/(m K)].
% [1] Cheong, J., Eun, Y., and Park, S.-Y. (2026).
% "Design of a Passive Sun-Pointing Mechanism for CubeSat Solar Panels."
% Aerospace, 13(7), 622. Table 4, FR4 row: 1900, 1150, 0.294.
% https://doi.org/10.3390/aerospace13070622
%
% [2] Park, Y.-K., Kim, G.-N., and Park, S.-Y. (2021).
% "Novel Structure and Thermal Design and Analysis for CubeSats in Formation Flying."
% Aerospace, 8(6), 150. Table 4, FR-4 row: 1850, 1300, 0.3.
% https://doi.org/10.3390/aerospace8060150
%
% [3] Kang, S.-J., and Oh, H.-U. (2016).
% "On-Orbit Thermal Design and Validation of 1 U Standardized CubeSat of STEP Cube Lab."
% International Journal of Aerospace Engineering, 2016, Article 4213189.
% Table 5, PCB/FR4 row: 1900, 1200, 0.1.
% https://doi.org/10.1155/2016/4213189
%
% [4] Mauro, S. (2016). "Thermal Analysis of Iodine Satellite (iSAT) from
% Preliminary Design Review (PDR) to Critical Design Review (CDR)."
% International Conference on Environmental Systems; NASA report M16-5317.
% NASA NTRS 20160009727, p. 17, Solar Panel FR4 row: 1850, 600, 0.3.
% https://ntrs.nasa.gov/citations/20160009727

%% Adopted property means
% Equal-weight arithmetic means calculated here from [1-4]:
% rho = (1900 + 1850 + 1900 + 1850)/4 = 1875;    adopted: 1900 kg/m^3.
% cp  = (1150 + 1300 + 1200 + 600)/4 = 1062.5;   adopted: 1100 J/(kg K).
% k   = (0.294 + 0.3 + 0.1 + 0.3)/4 = 0.2485;   adopted: 0.25 W/(m K).
% These are our averages of four published model inputs, not population means of flight arrays.
% Constant properties and a scalar through-thickness k are equivalent-plate approximations.

%% Carrier-thickness sources
% [5] Bhattarai, S., Kim, H., and Oh, H.-U. (2020). "CubeSat's Deployable Solar Panel with
% Viscoelastic Multilayered Stiffener for Launch Vibration Attenuation."
% International Journal of Aerospace Engineering, 2020, Article 8820619.
% Section 2.2, p. 3: 320 x 82 x 1.6 mm FR-4 carrier; rear stiffeners are additional.
% https://doi.org/10.1155/2020/8820619
%
% [6] Bhattarai, S., Go, J.-S., Kim, H., and Oh, H.-U. (2021).
% "Development of a Novel Deployable Solar Panel and Mechanism for
% 6U CubeSat of STEP Cube Lab-II." Aerospace, 8(3), 64.
% Section 3.1: 325.4 x 193 x 1.6 mm FR-4 carrier; rear stiffeners are additional.
% https://doi.org/10.3390/aerospace8030064
%
% Both carrier examples use 1.6 mm, so their two-example mean is 1.6 mm.
% This selected carrier thickness excludes cells, coverglass, adhesives and stiffeners.

%% Radiation constant and sink references
% [SIGMA] NIST, CODATA recommended values of the fundamental physical constants.
% Stefan-Boltzmann constant: 5.670374419...e-8 W/(m^2 K^4); code retains shown digits.
% https://physics.nist.gov/cuu/Constants/Table/allascii.txt
%
% [SPACE] NASA GSFC, Cosmicopia, "Heat, Temperature, and the Electromagnetic Spectrum."
% Entries "Temperature of the Universe" and "The Cold of Space" discuss the 2.7 K background.
% The adopted 3 K is a rounded ideal radiation sink; Earth and other local heating are omitted.
% https://cosmicopia.gsfc.nasa.gov/qa_sp_ht.html

%% Solar-irradiance reference
% [SUN] International Astronomical Union (2015), Resolution B3,
% "Recommended nominal conversion constants for selected solar and planetary properties."
% Nominal solar irradiance is 1361 W/m^2 at 1 AU, matching the supplied config.m value.
% https://www.iau.org/common/Uploaded%20files/IAUGA2015-Resolution-B3-recommended-nominal-conversion.pdf
%
% Existing solarFlux.m supplies S = cfg.S0*(cfg.AU/r_sun)^2*out.illum.
% The inverse-square factor follows conservation over spherical area; eclipse is applied upstream.
% [CELL]'s test irradiance normalizes gamma; runtime irradiance comes from solarFlux.m.

%% Incidence references
% [AOI] Sandia National Laboratories, PV Performance Modeling Collaborative, "POA Beam."
% Direct-beam projection: S_front = S_normal*cos(theta), with theta measured from the front normal
% toward the Sun; photon travel is opposite that Sun direction.
% https://pvpmc.sandia.gov/modeling-guide/1-weather-design-inputs/plane-of-array-poa-irradiance/calculating-poa-irradiance/poa-beam/
%
% MathWorks, "cosd": degree-based cosine evaluates cosd(90)=0 without a special power condition.
% https://www.mathworks.com/help/matlab/ref/double.cosd.html

%% Numerical-method references
% MathWorks MATLAB documentation: ode45 integrates the heat-balance ODE; odeset supplies
% tolerances and step limits, and deval evaluates the solution at the supplied orbit times.
% https://www.mathworks.com/help/matlab/ref/ode45.html
% https://www.mathworks.com/help/matlab/ref/odeset.html
% https://www.mathworks.com/help/matlab/ref/deval.html
%
% MathWorks, "interp1" and "trapz": linear irradiance interpolation and trapezoidal energy.
% Window-average power is integrated energy divided by elapsed time; 3600 converts J to Wh.
% https://www.mathworks.com/help/matlab/ref/double.interp1.html
% https://www.mathworks.com/help/matlab/ref/trapz.html

%% Graphing and input references
% MathWorks, "scatter": per-point color data map each angle/time sample to the angle color scale.
% https://www.mathworks.com/help/matlab/ref/scatter.html
% MathWorks, "turbo" and "flipud": reversing colormap rows changes the angle-to-color mapping.
% https://www.mathworks.com/help/matlab/ref/turbo.html
% https://www.mathworks.com/help/matlab/ref/flipud.html
% MathWorks, "nargin" and "validateattributes": optional-argument handling and input validation.
% https://www.mathworks.com/help/matlab/ref/nargin.html
% https://www.mathworks.com/help/matlab/ref/validateattributes.html

%% Project electrical assumptions
% Retained design estimates: cellEfficiency=0.30, cellCoverage=0.80, selfShadowRetained=1.00;
% mismatch=0.98, optical=0.97, wiring/diode=0.97, MPPT=0.98 and converter=0.95.
% These factors have no verified assembly-specific source or published average in this model.
%
% agingRetained=0.90 is a planning allowance without a mission duration or radiation fluence.
% It scales the BOL power result; it does not model a separately aged thermal state.

%% Project thermal and pointing assumptions
% Initial temperature 127 C is a chosen hot start; the front radiates with F_space=1 in sun/eclipses.
% Back/edges are adiabatic; Earth IR/albedo, mounting conduction and convection are omitted.
% MPPT shortfall remains panel heat; wiring/converter dissipation is assigned outside the plate.
%
% Default baseline pointing is normal to the Sun; the sweep defaults to A=0.20 m^2 and 0:1:90 deg.
% Direct sunlight is collimated; optical properties stay constant with incidence and edges absorb none.
% Default solver limits (10 s, RelTol=1e-7, AbsTol=1e-6 K) are numerical choices, not material data.
