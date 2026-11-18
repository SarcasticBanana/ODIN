function result = PMADThermalModel(inputs, overrides)
% PMADTHERMALMODEL Electrical conversion/distribution losses as heat [W].
% Use actual routed power, not the maximum available solar/battery power.
% No arguments evaluates the adjustable solar-powered load below.
% main_demo/componentThermalOrbit supplies separate nominal and peak histories, including actual charging.

%% Usage and input key
% m = finalPowerEPS.PMADThermalModel();  m.Qgen_W is heat inside the PMAD boundary.
% in = struct('time_s', [0 60 120], 'solarToLoad_W', [0 23.65 23.65], ...
%     'batteryToLoad_W', [23.65 0 0], 'batteryCharge_W', [0 10 10]);
% m = finalPowerEPS.PMADThermalModel(in);
% solarToLoad_W/batteryToLoad_W are delivered bus power; batteryCharge_W is
% electrical power delivered to battery terminals, NOT the stored-energy rate.
% Optional signals: controlPower_W, shuntPower_W and extraSwitchCurrent_A.
% Scalars expand in time; use time_s samples OR timeEdges_s for N interval values.
% Omitted signals use the top-of-file defaults; set inactive paths to zero.

% Shared defaults and assumptions: edit finalPowerEPS.modelSettings('pmad').
p=finalPowerEPS.modelSettings('pmad');

%% Parameter and input checks
if nargin < 1 || isempty(inputs), inputs = struct(); end
if nargin < 2, overrides = struct(); end
checkStruct(inputs, 'inputs');
checkStruct(overrides, 'overrides');
names = fieldnames(overrides);
for k = 1:numel(names)
    if ~isfield(p, names{k}), error('PMADThermalModel:parameter', 'Unknown parameter: %s', names{k}); end
    p.(names{k}) = overrides.(names{k});
end
efficiencies = {'mpptEfficiency','solarConverterEfficiency', ...
    'batteryDischargeEfficiency','batteryChargeConverterEfficiency'};
for k = 1:numel(efficiencies)
    v = p.(efficiencies{k});
    validateattributes(v, {'numeric'}, {'real','finite','scalar','>',0,'<=',1});
end
nonnegative = {'extraSwitchResistance_Ohm','defaultSolarToLoad_W','defaultBatteryToLoad_W', ...
    'defaultBatteryCharge_W','defaultControlPower_W','defaultShuntPower_W','balanceTolerance_W'};
for k = 1:numel(nonnegative)
    validateattributes(p.(nonnegative{k}), {'numeric'}, {'real','finite','scalar','nonnegative'});
end
validateattributes(p.defaultExtraSwitchCurrent_A, {'numeric'}, {'real','finite','scalar'});
p.showSummary = flagCheck(p.showSummary, 'showSummary');
p.warnAssumptions = flagCheck(p.warnAssumptions, 'warnAssumptions');
allowed = {'time_s','timeEdges_s','solarToLoad_W','batteryToLoad_W','batteryCharge_W', ...
    'controlPower_W','shuntPower_W','extraSwitchCurrent_A'};
names = fieldnames(inputs);
for k = 1:numel(names)
    if ~ismember(names{k}, allowed), error('PMADThermalModel:input', 'Unknown input: %s', names{k}); end
end
[t, edges, n] = timeBasis(inputs);
solar = signal(inputs, 'solarToLoad_W', p.defaultSolarToLoad_W, n);
discharge = signal(inputs, 'batteryToLoad_W', p.defaultBatteryToLoad_W, n);
charge = signal(inputs, 'batteryCharge_W', p.defaultBatteryCharge_W, n);
control = signal(inputs, 'controlPower_W', p.defaultControlPower_W, n);
shunt = signal(inputs, 'shuntPower_W', p.defaultShuntPower_W, n);
Iswitch = signal(inputs, 'extraSwitchCurrent_A', p.defaultExtraSwitchCurrent_A, n);
if any(solar < 0 | discharge < 0 | charge < 0 | control < 0 | shunt < 0)
    error('PMADThermalModel:power', 'All routed and dissipated powers must be nonnegative.');
end
if any(charge > 0 & discharge > 0)
    error('PMADThermalModel:batteryFlow', ...
        'A single modeled pack cannot have simultaneous net charge and discharge; supply actual net flows.');
end

%% PMAD heat generation section
% eta=Pout/Pin gives Qloss=Pout*(1/eta-1), not Pout*(1-eta) [TI].
% Charging includes the assumed MPPT path; the three routed output powers are disjoint.
etaSolar = p.mpptEfficiency * p.solarConverterEfficiency;
etaCharge = p.mpptEfficiency * p.batteryChargeConverterEfficiency;
Qsolar = solar .* (1 / etaSolar - 1);
Qdischarge = discharge .* (1 / p.batteryDischargeEfficiency - 1);
Qcharge = charge .* (1 / etaCharge - 1);
Qswitch = Iswitch.^2 .* p.extraSwitchResistance_Ohm;
Qgen = Qsolar + Qdischarge + Qcharge + Qswitch + control + shunt;

%% Electrical balance check
% Housekeeping, extra switch losses and shunt power require additional input power.
% Battery cell losses belong in BatteryThermalModel, not in these path efficiencies.
solarPathInput = solar ./ etaSolar + charge ./ etaCharge;
batteryPathInput = discharge ./ p.batteryDischargeEfficiency;
totalInput = solarPathInput + batteryPathInput + Qswitch + control + shunt;
totalOutput = solar + discharge + charge;
residual = totalInput - totalOutput - Qgen;
if any(~isfinite(Qgen) | ~isfinite(totalInput) | ~isfinite(totalOutput))
    error('PMADThermalModel:numerics', 'Calculated powers are not finite; check efficiencies and inputs.');
end
if any(abs(residual) > p.balanceTolerance_W + 32 * eps(max(totalInput, 1)))
    error('PMADThermalModel:balance', 'The PMAD electrical energy balance did not close.');
end
if p.warnAssumptions
    warning('PMADThermalModel:assumptions', ...
        ['Path efficiencies are provisional. Housekeeping/shunt/extra switches are included only as specified; ' ...
        'use actual routed power and avoid counting losses already included in an efficiency.']);
end

%% Results and command-window key
result.time_s = t;
result.timeEdges_s = edges;
result.Qgen_W = Qgen;
result.QsolarPath_W = Qsolar;
result.QbatteryDischargePath_W = Qdischarge;
result.QbatteryChargePath_W = Qcharge;
result.QextraSwitch_W = Qswitch;
result.Qcontrol_W = control;
result.Qshunt_W = shunt;
result.solarToLoad_W = solar;
result.batteryToLoad_W = discharge;
result.batteryCharge_W = charge;
result.solarPathInput_W = solarPathInput;
result.batteryPathInput_W = batteryPathInput;
result.additionalInput_W = Qswitch + control + shunt;
result.totalElectricalInput_W = totalInput;
result.totalElectricalOutput_W = totalOutput;
result.balanceResidual_W = residual;
result.solarPathEfficiency = etaSolar;
result.chargePathEfficiency = etaCharge;
result.summary = heatSummary(Qgen, t, edges);
result.params = p;
result.limitations = ['Assumed efficiencies until hardware curves are available; no controller, ' ...
    'battery SOC or component temperature solution. Power paths exclude battery internal heat ' ...
    'and solar-panel internal heat. Additional input power has no assigned source.'];
if p.showSummary
    if result.summary.duration_s == 0, fprintf('\nPMAD: single operating point; no orbit history supplied.\n'); end
    fprintf('\nPMAD thermal model: mean Qgen %.4g W; min %.4g W; max %.4g W.\n', ...
        result.summary.meanQgen_W, result.summary.minQgen_W, result.summary.maxQgen_W);
    fprintf('Key: Qgen=sum of conversion, extra-switch, control and shunt heat [W].\n');
    fprintf('Input powers describe actual routed electrical power; batteryCharge_W is at pack terminals.\n');
end
end

%% Local input and time helpers
function checkStruct(s, name)
if ~isstruct(s) || ~isscalar(s), error('PMADThermalModel:input', '%s must be a scalar struct.', name); end
end
function x = flagCheck(x, name)
if ~(islogical(x) || isnumeric(x)) || ~isscalar(x) || ~isreal(x) || ~isfinite(x) || ~(x == 0 || x == 1)
    error('PMADThermalModel:parameter', '%s must be true/false or 0/1.', name);
end
x = logical(x);
end
function [t, edges, n] = timeBasis(s)
edges = [];
if isfield(s, 'time_s') && isfield(s, 'timeEdges_s')
    error('PMADThermalModel:time', 'Supply time_s OR timeEdges_s.');
elseif isfield(s, 'timeEdges_s')
    edges = s.timeEdges_s;
    validateattributes(edges, {'numeric'}, {'real','finite','vector','nonempty'});
    edges = double(edges(:).');
    if numel(edges) < 2 || any(diff(edges) <= 0)
        error('PMADThermalModel:time', 'timeEdges_s needs at least two strictly increasing edges.');
    end
    t = edges(1:end-1) + diff(edges) / 2;
else
    if isfield(s, 'time_s'), t = s.time_s; else, t = 0; end
    validateattributes(t, {'numeric'}, {'real','finite','vector','nonempty'});
    t = double(t(:).');
    if any(diff(t) <= 0), error('PMADThermalModel:time', 'time_s must strictly increase.'); end
end
n = numel(t);
end
function x = signal(s, name, default, n)
if isfield(s, name), x = s.(name); else, x = default; end
validateattributes(x, {'numeric'}, {'real','finite','vector','nonempty'}, mfilename, name);
x = double(x(:).');
if isscalar(x), x = repmat(x, 1, n); end
if numel(x) ~= n, error('PMADThermalModel:length', '%s must be scalar or match the time basis.', name); end
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
% [TI] A. Castaldo, Switching Regulator Fundamentals, SNVA559C (2019), section 3.5.
% https://www.ti.com/lit/an/snva559c/snva559c.pdf
% Efficiency and conservation of electrical energy establish Pin=Pout+Qloss.
% MPPT 0.98 and converter 0.95 retain project assumptions; no PMAD datasheet is claimed.
% P=I^2*R applies only to explicitly separate switch/interconnect resistances.
