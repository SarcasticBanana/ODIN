function result = BatteryThermalModel(inputs, overrides)
% BATTERYTHERMALMODEL Cell/pack heat generation [W], not a temperature solver.
% Positive pack current/power means discharge; negative means charge.
% No arguments evaluates the adjustable representative operating point below.
% main_demo/componentThermalOrbit supplies separate nominal and peak orbit histories.

%% Usage and input key
% b = finalPowerEPS.BatteryThermalModel();  b.Qgen_W is total heat inside the modeled pack.
% in = struct('time_s', [0 60 120], 'packCurrent_A', [3 0 -2]);
% b = finalPowerEPS.BatteryThermalModel(in, struct('parallelCells', 2));
% Supply packCurrent_A OR packTerminalPower_W; cellVoltage_V is terminal voltage.
% Optional signals: cellTemperature_C, cellResistance_Ohm, entropyCoefficient_V_K,
% cellOCVoltage_V, heaterPower_W and electronicsPower_W; scalars expand in time.
% time_s holds samples; alternatively timeEdges_s holds N+1 edges for N interval values.
% Omitted signals use the top-of-file defaults; currents must be feasible actual operating values.
% Pack-terminal power is net external power; pack-local heater/electronics consumption stays inside this boundary.

% Shared defaults and assumptions: edit finalPowerEPS.modelSettings('battery').
p=finalPowerEPS.modelSettings('battery');

%% Parameter and input checks
if nargin < 1 || isempty(inputs), inputs = struct(); end
if nargin < 2, overrides = struct(); end
checkStruct(inputs, 'inputs');
checkStruct(overrides, 'overrides');
p = applyOverrides(p, overrides);
positive = {'cellNominalCapacity_Ah','cellMinimumCapacity_Ah','cellNominalVoltage_V', ...
    'cellStandardChargeCurrent_A','cellChargeVoltage_V','cellChargeTerminationCurrent_A', ...
    'cellMaxChargeCurrent_A','cellStandardDischargeCurrent_A','cellDischargeCutoffVoltage_V', ...
    'cellMaxDischargeCurrent_A','cellMass_kg','cellDiameter_m','cellHeight_m', ...
    'sheetACImpedanceMax_Ohm'};
for k = 1:numel(positive), scalarCheck(p.(positive{k}), positive{k}, 0, Inf, true); end
scalarCheck(p.seriesCells, 'seriesCells', 1, Inf, false);
scalarCheck(p.parallelCells, 'parallelCells', 1, Inf, false);
if p.seriesCells ~= fix(p.seriesCells) || p.parallelCells ~= fix(p.parallelCells)
    error('BatteryThermalModel:layout', 'Cell counts must be positive integers.');
end
scalarCheck(p.defaultBusLoad_W, 'defaultBusLoad_W', 0, Inf, false);
scalarCheck(p.defaultBatteryToBusEfficiency, 'defaultBatteryToBusEfficiency', 0, 1, true);
scalarCheck(p.cellDCResistance_Ohm, 'cellDCResistance_Ohm', 0, Inf, false);
scalarCheck(p.packInterconnectResistance_Ohm, 'packInterconnectResistance_Ohm', 0, Inf, false);
scalarCheck(p.defaultHeaterPower_W, 'defaultHeaterPower_W', 0, Inf, false);
scalarCheck(p.defaultElectronicsPower_W, 'defaultElectronicsPower_W', 0, Inf, false);
scalarCheck(p.limitTolerance, 'limitTolerance', 0, Inf, false);
scalarCheck(p.defaultCellTemperature_C, 'defaultCellTemperature_C', -273.15, Inf, true);
scalarCheck(p.entropyCoefficient_V_K, 'entropyCoefficient_V_K', -Inf, Inf, false);
flags = {'cellDCResistanceIsMeasured','showSummary','warnAssumptions'};
for k = 1:numel(flags), p.(flags{k}) = flagCheck(p.(flags{k}), flags{k}); end
ranges = {'cellChargeTemperatureRange_C','cellDischargeTemperatureRange_C'};
for k = 1:numel(ranges)
    v = p.(ranges{k});
    validateattributes(v, {'numeric'}, {'real','finite','vector','numel',2});
    if v(1) <= -273.15 || v(2) < v(1)
        error('BatteryThermalModel:range', '%s must be ordered above absolute zero.', ranges{k});
    end
end
if p.cellMinimumCapacity_Ah > p.cellNominalCapacity_Ah || ...
        p.cellDischargeCutoffVoltage_V >= p.cellChargeVoltage_V
    error('BatteryThermalModel:specification', 'Check the capacity and voltage limits.');
end
allowed = {'time_s','timeEdges_s','packCurrent_A','packTerminalPower_W','cellVoltage_V', ...
    'cellTemperature_C','cellResistance_Ohm','entropyCoefficient_V_K', ...
    'cellOCVoltage_V','heaterPower_W','electronicsPower_W'};
checkFields(inputs, allowed);
[t, edges, n] = timeBasis(inputs);
V = signal(inputs, 'cellVoltage_V', p.cellNominalVoltage_V, n);
T_C = signal(inputs, 'cellTemperature_C', p.defaultCellTemperature_C, n);
R = signal(inputs, 'cellResistance_Ohm', p.cellDCResistance_Ohm, n);
entropy = signal(inputs, 'entropyCoefficient_V_K', p.entropyCoefficient_V_K, n);
heater = signal(inputs, 'heaterPower_W', p.defaultHeaterPower_W, n);
electronics = signal(inputs, 'electronicsPower_W', p.defaultElectronicsPower_W, n);
if any(V <= 0) || any(T_C <= -273.15) || any(R < 0) || any(heater < 0) || any(electronics < 0)
    error('BatteryThermalModel:physics', 'Require V>0, T>absolute zero, and nonnegative R/auxiliary power.');
end
if isfield(inputs, 'packCurrent_A') && isfield(inputs, 'packTerminalPower_W')
    error('BatteryThermalModel:current', 'Supply packCurrent_A OR packTerminalPower_W, not both.');
elseif isfield(inputs, 'packCurrent_A')
    I = signal(inputs, 'packCurrent_A', 0, n);
    currentSource = 'supplied signed pack current';
elseif isfield(inputs, 'packTerminalPower_W')
    P = signal(inputs, 'packTerminalPower_W', 0, n);
    I = currentForPower(P + heater + electronics, p.seriesCells .* V, p.packInterconnectResistance_Ohm);
    currentSource = 'supplied external terminal power plus pack auxiliaries and interconnect loss';
else
    P = repmat(p.defaultBusLoad_W / p.defaultBatteryToBusEfficiency,1,n);
    I = currentForPower(P + heater + electronics, p.seriesCells .* V, p.packInterconnectResistance_Ohm);
    currentSource = 'representative load and assumed battery-to-bus efficiency';
end

%% Battery heat generation section
% Equal current sharing is assumed; use measured DC resistance or OCV/terminal voltages when available.
% With discharge-positive current, reversible heat is -I*T*d(OCV)/dT [HEAT].
cellCount = p.seriesCells * p.parallelCells;
Icell = I ./ p.parallelCells;
Qirr = cellCount .* Icell.^2 .* R;
method = 'I^2 R, effective cell DC resistance';
if isfield(inputs, 'cellOCVoltage_V')
    if ~isfield(inputs, 'cellVoltage_V')
        error('BatteryThermalModel:voltage', 'OCV-based heat also requires cellVoltage_V.');
    end
    OCV = signal(inputs, 'cellOCVoltage_V', 0, n);
    Qirr = cellCount .* Icell .* (OCV - V);
    if any(OCV <= 0) || any(Qirr < -p.limitTolerance)
        error('BatteryThermalModel:voltage', 'OCV/terminal voltage and current imply negative irreversible heat.');
    end
    Qirr = max(Qirr, 0); % remove roundoff only; inconsistent voltages were rejected above
    method = 'I*(OCV-terminal voltage), no additional I^2 R term';
end
Qrev = -cellCount .* Icell .* (T_C + 273.15) .* entropy;
Qinterconnect = I.^2 .* p.packInterconnectResistance_Ohm;
Qgen = Qirr + Qrev + Qinterconnect + heater + electronics;
if any(~isfinite(Qgen))
    error('BatteryThermalModel:numerics', 'Calculated heat is not finite; check input magnitudes.');
end

%% Operating-limit checks
% These flags report excursions; this heat-source function does not implement a BMS or CC/CV charger.
tol = p.limitTolerance;
discharge = Icell > 0;
charge = Icell < 0;
currentExceeded = (discharge & Icell > p.cellMaxDischargeCurrent_A + tol) | ...
    (charge & -Icell > p.cellMaxChargeCurrent_A + tol);
temperatureExceeded = (charge & (T_C < p.cellChargeTemperatureRange_C(1) | ...
    T_C > p.cellChargeTemperatureRange_C(2))) | ...
    (discharge & (T_C < p.cellDischargeTemperatureRange_C(1) | ...
    T_C > p.cellDischargeTemperatureRange_C(2)));
voltageExceeded = V < p.cellDischargeCutoffVoltage_V - tol | V > p.cellChargeVoltage_V + tol;
provisionalR = ~isfield(inputs, 'cellOCVoltage_V') && ...
    ~isfield(inputs, 'cellResistance_Ohm') && ~p.cellDCResistanceIsMeasured;
if p.warnAssumptions && provisionalR
    warning('BatteryThermalModel:DCResistance', ...
        'Heat uses assumed DC resistance %.4g ohm/cell; the sheet specifies AC impedance, not DC resistance.', ...
        p.cellDCResistance_Ohm);
end
if p.warnAssumptions && any(currentExceeded | temperatureExceeded | voltageExceeded)
    warning('BatteryThermalModel:operatingLimits', ...
        'Some supplied operating points exceed cell current, voltage or temperature limits; inspect result flags.');
end

%% Results and command-window key
result.time_s = t;
result.timeEdges_s = edges;
result.Qgen_W = Qgen;                    % total pack heat source; signed if entropy is supplied
result.Qirreversible_W = Qirr;
result.Qreversible_W = Qrev;
result.Qinterconnect_W = Qinterconnect;
result.Qheater_W = heater;
result.Qelectronics_W = electronics;
result.cellQgen_W = (Qirr + Qrev) ./ cellCount;
result.packCurrent_A = I;
result.cellCurrent_A = Icell;
result.packTerminalVoltage_V = p.seriesCells .* V - I .* p.packInterconnectResistance_Ohm;
if any(result.packTerminalVoltage_V <= 0)
    error('BatteryThermalModel:terminalVoltage', 'Interconnect voltage drop leaves nonpositive pack-terminal voltage.');
end
result.packTerminalPower_W = I .* result.packTerminalVoltage_V - heater - electronics;
result.cellVoltage_V = V;
result.cellTemperature_C = T_C;
result.irreversibleMethod = method;
result.currentSource = currentSource;
result.usesProvisionalDCResistance = provisionalR;
result.currentLimitExceeded = currentExceeded;
result.temperatureLimitExceeded = temperatureExceeded;
result.voltageLimitExceeded = voltageExceeded;
result.withinOperatingLimits = ~(currentExceeded | temperatureExceeded | voltageExceeded);
result.pack.nominalVoltage_V = p.seriesCells * p.cellNominalVoltage_V;
result.pack.nominalCapacity_Ah = p.parallelCells * p.cellNominalCapacity_Ah;
result.pack.nominalEnergy_Wh = result.pack.nominalVoltage_V * result.pack.nominalCapacity_Ah; % V*Ah estimate
result.pack.minimumCapacityEnergy_Wh = p.seriesCells * p.cellNominalVoltage_V * ...
    p.parallelCells * p.cellMinimumCapacity_Ah;
result.pack.cellOnlyMass_kg = cellCount * p.cellMass_kg;
result.pack.maxCellLimitedDischarge_A = p.parallelCells * p.cellMaxDischargeCurrent_A;
result.pack.maxCellLimitedCharge_A = p.parallelCells * p.cellMaxChargeCurrent_A;
result.pack.effectiveDCResistance_Ohm = (p.seriesCells / p.parallelCells) .* R;
result.summary = heatSummary(Qgen, t, edges);
result.params = p;
result.limitations = ['Preliminary normal-operation heat source; no temperature, SOC, capacity fade, ' ...
    'cell imbalance or abuse model. Entropy defaults to zero; DC resistance is unvalidated. ' ...
    'Pack mass/limits exclude enclosure, BMS and interconnects.'];
if p.showSummary
    if result.summary.duration_s == 0, fprintf('\nBattery: single operating point; no orbit history supplied.\n'); end
    fprintf('\nBattery thermal model: mean Qgen %.4g W; min %.4g W; max %.4g W.\n', ...
        result.summary.meanQgen_W, result.summary.minQgen_W, result.summary.maxQgen_W);
    fprintf('Pack assumption: %gS%gP, %.3g V, %.3g Ah, %.4g Wh nominal; %.4g kg cells only.\n', ...
        p.seriesCells, p.parallelCells, result.pack.nominalVoltage_V, ...
        result.pack.nominalCapacity_Ah, result.pack.nominalEnergy_Wh, result.pack.cellOnlyMass_kg);
    fprintf('Key: Qgen=Qirreversible+Qreversible+Qinterconnect+Qheater+Qelectronics [W].\n');
    fprintf('Positive current/power means discharge; negative means charge.\n');
    fprintf('Current source: %s.\n', currentSource);
    fprintf('Entropy heat is omitted when its coefficient is zero; no temperature is solved.\n');
end
end

%% Local input helpers
function I = currentForPower(P, Vstring, R)
% Solve P=I*Vstring-I^2*R on the stable branch; negative P naturally gives charging current.
D = Vstring.^2 - 4 .* R .* P;
if any(D < 0), error('BatteryThermalModel:power', 'Requested power exceeds the interconnect-resistance transfer limit.'); end
I = 2 .* P ./ (Vstring + sqrt(D));
end
function checkStruct(s, name)
if ~isstruct(s) || ~isscalar(s), error('BatteryThermalModel:input', '%s must be a scalar struct.', name); end
end
function p = applyOverrides(p, changes)
names = fieldnames(changes);
for k = 1:numel(names)
    if ~isfield(p, names{k}), error('BatteryThermalModel:parameter', 'Unknown parameter: %s', names{k}); end
    p.(names{k}) = changes.(names{k});
end
end
function scalarCheck(x, name, lo, hi, strictLo)
validateattributes(x, {'numeric'}, {'real','finite','scalar'}, mfilename, name);
if x < lo || x > hi || (strictLo && x == lo)
    error('BatteryThermalModel:parameter', 'Invalid parameter: %s', name);
end
end
function x = flagCheck(x, name)
if ~(islogical(x) || isnumeric(x)) || ~isscalar(x) || ~isreal(x) || ~isfinite(x) || ~(x == 0 || x == 1)
    error('BatteryThermalModel:parameter', '%s must be true/false or 0/1.', name);
end
x = logical(x);
end
function checkFields(s, allowed)
names = fieldnames(s);
for k = 1:numel(names)
    if ~ismember(names{k}, allowed), error('BatteryThermalModel:input', 'Unknown input: %s', names{k}); end
end
end
function [t, edges, n] = timeBasis(s)
edges = [];
if isfield(s, 'time_s') && isfield(s, 'timeEdges_s')
    error('BatteryThermalModel:time', 'Supply time_s OR timeEdges_s.');
elseif isfield(s, 'timeEdges_s')
    edges = s.timeEdges_s;
    validateattributes(edges, {'numeric'}, {'real','finite','vector','nonempty'});
    edges = double(edges(:).');
    if numel(edges) < 2 || any(diff(edges) <= 0)
        error('BatteryThermalModel:time', 'timeEdges_s needs at least two strictly increasing edges.');
    end
    t = edges(1:end-1) + diff(edges) / 2;
else
    if isfield(s, 'time_s'), t = s.time_s; else, t = 0; end
    validateattributes(t, {'numeric'}, {'real','finite','vector','nonempty'});
    t = double(t(:).');
    if any(diff(t) <= 0), error('BatteryThermalModel:time', 'time_s must strictly increase.'); end
end
n = numel(t);
end
function x = signal(s, name, default, n)
if isfield(s, name), x = s.(name); else, x = default; end
validateattributes(x, {'numeric'}, {'real','finite','vector','nonempty'}, mfilename, name);
x = double(x(:).');
if isscalar(x), x = repmat(x, 1, n); end
if numel(x) ~= n, error('BatteryThermalModel:length', '%s must be scalar or match the time basis.', name); end
end
function s = heatSummary(q, t, edges)
s.minQgen_W = min(q); s.maxQgen_W = max(q);
if ~isempty(edges)
    s.duration_s = edges(end) - edges(1);
    s.heatEnergy_J = sum(q .* diff(edges));
elseif numel(t) > 1
    s.duration_s = t(end) - t(1);
    s.heatEnergy_J = trapz(t, q);
else
    s.duration_s = 0; s.heatEnergy_J = NaN;
end
if s.duration_s > 0, s.meanQgen_W = s.heatEnergy_J / s.duration_s; else, s.meanQgen_W = q(1); end
s.heatEnergy_Wh = s.heatEnergy_J / 3600;
end

%% References
% [BAT] AA Portable Power Corp., LC-18650-2600, VER A, 2005-01-03.
% https://www.batteryspace.com/prod-specs/18650_2600.pdf
% [NASA] Cooper and London, Practical Battery Thermal Modeling Techniques,
% TFAWS 2020, slide 4: reduced I^2 R approximation and its limitations.
% https://ntrs.nasa.gov/citations/20205006303
% [HEAT] Paccha-Herrera et al., Batteries 6(3), 40 (2020), section 2.1.1.
% https://doi.org/10.3390/batteries6030040
% For discharge-positive current: Qirr=I*(OCV-V), Qrev=-I*T*d(OCV)/dT.
% [AC/DC] Li-Ion BMS, Resistance vs impedance: 1 kHz impedance is not DC resistance.
% https://liionbms.com/php/wp_resistance_vs_impedance.php
