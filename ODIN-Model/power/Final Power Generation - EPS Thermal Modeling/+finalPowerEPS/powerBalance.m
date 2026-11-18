function result = powerBalance(cfg, overrides)
% POWERBALANCE Worst-shadow solar/battery supply against a continuous all-on load.
% Run result = finalPowerEPS.powerBalance(), or let main_demo.m call it after its existing plots.
% Required files: config.m, solarFlux.m and the updated arrayPower.m.

% Shared defaults and assumptions: edit finalPowerEPS.modelSettings('powerBalance').
s=finalPowerEPS.modelSettings('powerBalance');

%% Parameter checks
if nargin < 1 || isempty(cfg)
    cfg = finalPowerEPS.config();
end
if nargin < 2
    overrides = struct();
end
validateattributes(cfg, {'struct'}, {'scalar'});
validateattributes(overrides, {'struct'}, {'scalar'});
names = fieldnames(overrides);
for k = 1:numel(names)
    if ~isfield(s, names{k})
        error('powerBalance:parameter', 'Unknown setting: %s', names{k});
    end
    s.(names{k}) = overrides.(names{k});
end
validateattributes(s.arrayOverrides, {'struct'}, {'scalar'});
scalarFields = setdiff(fieldnames(s), {'arrayOverrides','arrayArea_m2', ...
    'batteryMaxCharge_W','figureColor','figurePosition','eclipseShadeColor', ...
    'batteryTradeFigureSize_px','batteryOperationFigureSize_px', ...
    'balanceColors','balanceLineStyles','balanceLabels','batteryTradeDischarge_percent', ...
    'batteryTradeColorAnchors','batteryTradeView_deg','batteryTradeMeshColor','batteryTradeCurveCommands_percent'});
for k = 1:numel(scalarFields)
    validateattributes(s.(scalarFields{k}), {'numeric','logical'}, ...
        {'real','finite','scalar'});
end
flagFields = {'allowBatteryInSunlight','showSummary','makePlot','warnModelLimits','runBatteryTrade','makeBatteryTradePlots','makeBatteryTradeSurfaceViews'};
for k = 1:numel(flagFields)
    name = flagFields{k};
    validateattributes(double(s.(name)), {'numeric'}, {'integer','>=',0,'<=',1});
    s.(name) = logical(s.(name));
end
positiveFields = {'sunDistance_AU','continuousLoad_W','outputStep_s', ...
    'batteryCapacityBOL_Wh','batteryCapacityRetained','batteryMaxDischarge_W', ...
    'batteryToBusEfficiency','thermalPeriodicTolerance_C','lineWidth','balanceTickStep_percent', ...
    'batteryTradeLineWidth','batteryTradeHighlightWidth', ...
    'batteryTradeMeshRateStep_percent','batteryTradeMeshTimeStep_min','batteryTradeMeshLineWidth'};
for k = 1:numel(positiveFields)
    validateattributes(s.(positiveFields{k}), {'numeric'}, {'positive'});
end
fractionFields = {'batteryCapacityRetained','batteryInitialSOC','batteryMinimumSOC', ...
    'batteryMaximumSOC','batteryToBusEfficiency','batteryChargeEfficiency','eclipseShadeAlpha','batteryTradeSurfaceAlpha'};
for k = 1:numel(fractionFields)
    validateattributes(s.(fractionFields{k}), {'numeric'}, {'>=',0,'<=',1});
end
validateattributes(s.balanceColors, {'numeric'}, {'real','finite','size',[5 3],'>=',0,'<=',1});
validateattributes(s.balanceAxisPadding, {'numeric'}, {'>=',1});
validateattributes(s.supplementalMargin_percent, {'numeric'}, {'nonnegative'});
validateattributes(s.balanceLabels, {'cell'}, {'vector','numel',5});
validateattributes(s.balanceLineStyles, {'cell'}, {'vector','numel',5});
for k = 1:5
    validateattributes(s.balanceLabels{k}, {'char'}, {'row','nonempty'});
    s.balanceLineStyles{k} = validatestring(s.balanceLineStyles{k}, {'-','--',':','-.'});
end
validateattributes(s.batteryTradeDischarge_percent, {'numeric'}, ...
    {'real','finite','vector','nonempty','>=',0,'<=',100});
s.batteryTradeDischarge_percent = unique(s.batteryTradeDischarge_percent(:));
validateattributes(s.batteryTradeColorLevels, {'numeric'}, {'integer','>=',2});
validateattributes(s.batteryTradeColorAnchors, {'numeric'}, ...
    {'real','finite','2d','ncols',3,'>=',0,'<=',1});
if size(s.batteryTradeColorAnchors,1) < 2
    error('powerBalance:tradeColors', 'Need at least two color anchors.');
end
validateattributes(s.batteryTradeCurveCommands_percent,{'numeric'},{'real','finite','vector','nonempty','>=',0,'<=',100});
validateattributes(s.batteryTradeView_deg, {'numeric'}, {'real','finite','vector','numel',2});
validateattributes(s.batteryTradeTolerance_percent, {'numeric'}, {'nonnegative'});
validateattributes(s.orbitAltitude_km, {'numeric'}, {'nonnegative'});
validateattributes(s.eclipseDuration_min, {'numeric'}, {'positive'});
validateattributes(s.thermalWarmupMaxOrbits, {'numeric'}, {'integer','positive'});
validateattributes(s.energyTolerance_Wh, {'numeric'}, {'nonnegative'});
validateattributes(s.batteryMaxCharge_W, {'numeric'}, {'real','scalar','nonnegative','nonnan'});
if s.batteryMinimumSOC >= s.batteryMaximumSOC || ...
        s.batteryInitialSOC < s.batteryMinimumSOC || s.batteryInitialSOC > s.batteryMaximumSOC
    error('powerBalance:SOC', 'Need minimum SOC <= initial SOC <= maximum SOC, with min < max.');
end
for name = {'Re','mu','AU','S0'}
    if ~isfield(cfg, name{1})
        error('powerBalance:config', 'cfg must contain %s.', name{1});
    end
    validateattributes(cfg.(name{1}), {'numeric'}, {'real','finite','scalar','positive'});
end
for name = {'figureColor','eclipseShadeColor','batteryTradeMeshColor'}
    validateattributes(s.(name{1}), {'numeric'}, {'real','finite','vector','numel',3,'>=',0,'<=',1});
    color = s.(name{1});
    s.(name{1}) = color(:).';
end
validateattributes(s.figurePosition, {'numeric'}, {'real','finite','vector','numel',4});
s.figurePosition = s.figurePosition(:).';
validateattributes(s.figurePosition(3:4), {'numeric'}, {'positive'});
for name={'batteryTradeFigureSize_px','batteryOperationFigureSize_px'}
    validateattributes(s.(name{1}),{'numeric'},{'real','finite','positive','vector','numel',2});
end

%% Array settings section
% Keep area in one place: arrayPower's selectedArrayArea_m2 unless explicitly overridden here.
arrayDefaults = finalPowerEPS.arrayPower();
arraySettings = s.arrayOverrides;
if isempty(s.arrayArea_m2)
    s.arrayArea_m2 = arrayDefaults.selectedArrayArea_m2;
    if isfield(arraySettings, 'selectedArrayArea_m2')
        s.arrayArea_m2 = arraySettings.selectedArrayArea_m2;
    end
end
validateattributes(s.arrayArea_m2, {'numeric'}, {'real','finite','scalar','positive'});
if ~isfield(arraySettings, 'initialTemperature_C')
    arraySettings.initialTemperature_C = arrayDefaults.initialTemperature_C;
end
arraySettings.runIncidenceSweep = false;
arraySettings.makePlots = false;
arraySettings.makeMainArrayPlots = false;
arraySettings.showAreaSummary = false;
arraySettings.showIncidenceSummary = false;
arraySettings.warnTemperatureExtrapolation = false;
arraySettings.warnLumpedModel = false;
arraySettings.warnProductTemperature = false;
arraySettings.irradianceInterpolation = 'previous'; % abrupt envelope boundary, no linear twilight

if ~isfield(arraySettings,'initialBackTemperature_C')
    arraySettings.initialBackTemperature_C=arrayDefaults.initialBackTemperature_C;
end

%% Worst-shadow orbit section
% Eclipse begins at t=0 so a full battery can discharge, then recharge in the following sunlight.
% The synthetic illumination history uses solarFlux; it does not replace the team's propagator.
period_s = 2*pi*sqrt((cfg.Re + s.orbitAltitude_km)^3 / cfg.mu);
eclipse_s = s.eclipseDuration_min * 60;
if eclipse_s >= period_s
    error('powerBalance:orbit', 'Eclipse duration must be shorter than the orbital period.');
end
t = unique([0:s.outputStep_s:period_s, eclipse_s, period_s]);
envelope.t = t;
envelope.illum = double(t >= eclipse_s);
envelope.rsun = repmat([cfg.AU*s.sunDistance_AU; 0; 0], 1, numel(t));
envelope.T = period_s;
% Environmental geometry for the independent envelope: beta-zero circular orbit,
% centered on its prescribed eclipse. Illumination still follows the original
% 35.8-minute design envelope, not a propagated event/date.
phase=pi+2*pi*(t-eclipse_s/2)/period_s;
envelope.r=(cfg.Re+s.orbitAltitude_km)*[cos(phase);sin(phase);zeros(size(t))];
envelope.v=(cfg.Re+s.orbitAltitude_km)*2*pi/period_s*[-sin(phase);cos(phase);zeros(size(t))];


%% Transient heat section
% Reuse the coupled electrical/thermal equations, including their existing loss assumptions.
% Warm-up affects only this independent worst-shadow study.
periodicError_C = Inf;
for orbit = 1:s.thermalWarmupMaxOrbits
    array = finalPowerEPS.arrayPower(envelope, cfg, s.arrayArea_m2, arraySettings);
    periodicError_C = max(abs([array.cellTemperature_C(end)-arraySettings.initialTemperature_C, ...
        array.backTemperature_C(end)-arraySettings.initialBackTemperature_C]));
    if periodicError_C <= s.thermalPeriodicTolerance_C
        break
    end
    arraySettings.initialTemperature_C = array.cellTemperature_C(end);
    arraySettings.initialBackTemperature_C = array.backTemperature_C(end);
end
if s.warnModelLimits && periodicError_C > s.thermalPeriodicTolerance_C
    warning('powerBalance:thermalSettling', 'Thermal orbit did not settle within %g C in %d repeats.', ...
        s.thermalPeriodicTolerance_C, s.thermalWarmupMaxOrbits);
end

%% Power input section
% Average each smooth sunlit interval; the shadow boundary is an exact time-grid point.
% Planning bus power retains arrayPower's assumed delivery factors and 10% aging allowance.
dt = diff(t);
sunlit = envelope.illum(1:end-1) > 0;
solarPoint_W = array.busPlanning_W(1,:);
solar_W = (solarPoint_W(1:end-1) + solarPoint_W(2:end))/2;
solar_W(~sunlit) = 0;
load_W = s.continuousLoad_W * ones(size(dt));
deficit_W = max(load_W - solar_W, 0);
surplus_W = max(solar_W - load_W, 0);
batteryAllowed = ~sunlit | s.allowBatteryInSunlight;
requestedBus_W = deficit_W .* batteryAllowed;

%% Battery energy section
% Actual battery output follows demand and is limited by both power and stored energy.
% Charging uses only solar surplus; the floor stops discharge instead of hiding a load shortfall.
capacity_Wh = s.batteryCapacityBOL_Wh * s.batteryCapacityRetained;
floor_Wh = capacity_Wh * s.batteryMinimumSOC;
ceiling_Wh = capacity_Wh * s.batteryMaximumSOC;
stored_Wh = zeros(size(t));
stored_Wh(1) = capacity_Wh * s.batteryInitialSOC;
batteryBus_W = zeros(size(dt));
batteryPack_W = zeros(size(dt));
chargePack_W = zeros(size(dt));
potentialChargePack_W = min(surplus_W * s.batteryChargeEfficiency, s.batteryMaxCharge_W);
for k = 1:numel(dt)
    energyLimitedPack_W = max(stored_Wh(k) - floor_Wh, 0) * 3600 / dt(k);
    batteryPack_W(k) = min([requestedBus_W(k) / s.batteryToBusEfficiency, ...
        s.batteryMaxDischarge_W, energyLimitedPack_W]);
    batteryBus_W(k) = batteryPack_W(k) * s.batteryToBusEfficiency;
    roomLimitedPack_W = max(ceiling_Wh - stored_Wh(k), 0) * 3600 / dt(k);
    chargePack_W(k) = min(potentialChargePack_W(k), roomLimitedPack_W);
    stored_Wh(k+1) = stored_Wh(k) + (chargePack_W(k) - batteryPack_W(k)) * dt(k)/3600;
end
unserved_W = max(load_W - solar_W - batteryBus_W, 0);

%% Energy balance section
% Battery discharge is stored solar energy, not an additional source in the orbit energy balance.
% Recharge margin uses potential surplus before curtailment when a full battery cannot accept more.
integrateWh = @(power_W) sum(power_W .* dt)/3600;
solar_Wh = integrateWh(solar_W);
load_Wh = integrateWh(load_W);
requiredPackDischarge_Wh = integrateWh(requestedBus_W / s.batteryToBusEfficiency);
potentialPackCharge_Wh = integrateWh(potentialChargePack_W);
rechargeMargin_Wh = potentialPackCharge_Wh - requiredPackDischarge_Wh;
unserved_Wh = integrateWh(unserved_W);
uncoveredSunlight_Wh = integrateWh(deficit_W .* ~batteryAllowed);
energyPositive = rechargeMargin_Wh >= -s.energyTolerance_Wh && ...
    uncoveredSunlight_Wh <= s.energyTolerance_Wh;
loadCovered = unserved_Wh <= s.energyTolerance_Wh;
if energyPositive && loadCovered
    status = 'ENERGY POSITIVE under the stated assumptions';
elseif energyPositive
    status = 'ENERGY POSITIVE, but battery power/reserve cannot cover every load';
else
    status = 'ENERGY NEGATIVE under the stated assumptions';
end
surplus_Wh = integrateWh(surplus_W);
minimumChargeEfficiency = Inf;
if surplus_Wh > 0
    minimumChargeEfficiency = requiredPackDischarge_Wh / surplus_Wh;
end

%% Results section
% Power histories are interval averages (N-1); stored energy/SOC use the N time boundaries.
result.params = s;
result.arrayParams = array.params;
result.time_s = t;
result.intervalDuration_s = dt;
result.sunlit = sunlit;
result.solarAvailable_W = solar_W;
result.batteryToBus_W = batteryBus_W;
result.combinedAvailable_W = solar_W + batteryBus_W;
result.load_W = load_W;
result.unserved_W = unserved_W;
result.batteryPackDischarge_W = batteryPack_W;
result.batteryPackCharge_W = chargePack_W;
result.storedEnergy_Wh = stored_Wh;
result.batterySOC = stored_Wh / capacity_Wh;
result.cellTemperature_C = array.cellTemperature_C;
result.thermalPeriodicError_C = periodicError_C;
result.thermalWarmupOrbits = orbit;
result.energy.solar_Wh = solar_Wh;
result.energy.load_Wh = load_Wh;
result.energy.solarMinusLoad_Wh = solar_Wh - load_Wh;
result.energy.requiredPackDischarge_Wh = requiredPackDischarge_Wh;
result.energy.actualPackDischarge_Wh = integrateWh(batteryPack_W);
result.energy.actualPackCharge_Wh = integrateWh(chargePack_W);
result.energy.potentialPackCharge_Wh = potentialPackCharge_Wh;
result.energy.rechargeMargin_Wh = rechargeMargin_Wh;
result.energy.unserved_Wh = unserved_Wh;
result.energy.packToBusLoss_Wh = integrateWh(batteryPack_W - batteryBus_W);
result.energy.batteryStoredEnergyChange_Wh = stored_Wh(end) - stored_Wh(1);
result.energy.minimumChargeEfficiencyWithoutRateLimit = minimumChargeEfficiency;
result.energy.chargeRateLimitIncluded_W = s.batteryMaxCharge_W;
result.energyPositive = energyPositive;
result.loadCovered = loadCovered;
result.status = status;
result.orbitPeriod_min = period_s/60;
result.unshadowedIrradiance_W_m2 = cfg.S0 / s.sunDistance_AU^2;
result.batteryPowerLimitAtBus_W = s.batteryMaxDischarge_W * s.batteryToBusEfficiency;
result.summary = table(period_s/60, s.eclipseDuration_min, s.arrayArea_m2, ...
    s.continuousLoad_W, solar_Wh, load_Wh, rechargeMargin_Wh, ...
    min(result.batterySOC)*100, unserved_Wh, 'VariableNames', ...
    {'Orbit_min','Shadow_min','Area_m2','Load_W','Solar_Wh','Load_Wh', ...
     'RechargeMargin_Wh','MinimumSOC_percent','Unserved_Wh'});
result.summary.Properties.VariableDescriptions = { ...
    'Duration of the constructed worst-shadow orbit [min].'; ...
    'Duration treated as fully eclipsed [min], including the assumed penumbra allowance.'; ...
    'Total physical solar-array front-face area [m^2].'; ...
    'Continuous electrical demand at the load bus [W].'; ...
    'Available planning solar energy at the bus over this orbit, before charging or curtailment [Wh].'; ...
    'Electrical energy demanded by the load over this orbit [Wh].'; ...
    'Potential pack charging energy minus required pack discharge energy in the baseline dispatch [Wh].'; ...
    'Minimum battery state of charge in the baseline dispatch [% of retained capacity], not the comparison cases.'; ...
    'Load energy left unmet in the baseline dispatch because available supply is insufficient [Wh].'};

%% Solar power statistics section
% Orbit average includes eclipse; instantaneous minima are reported separately.
result.solarPowerStats.orbitAveragePlanning_W = solar_Wh * 3600 / period_s;
result.solarPowerStats.minimumPlanning_W = min(solar_W);
result.solarPowerStats.minimumSunlitPlanning_W = NaN;
if any(sunlit)
    result.solarPowerStats.minimumSunlitPlanning_W = min(solar_W(sunlit));
end

%% Comparison scenarios section
% Each case starts with its own identical battery state; positive margin is available headroom.
% The maximum-discharge case exports all battery output, with no simultaneous recharge.
zero_W = zeros(size(dt));
batteryOnly = comparisonBattery(load_W, zero_W, dt, s);
maximum = comparisonBattery(s.batteryMaxDischarge_W * s.batteryToBusEfficiency * ...
    ones(size(dt)), zero_W, dt, s);
conditional = comparisonBattery(load_W .* ~sunlit, surplus_W .* sunlit, dt, s);

% Supplemental discharge targets the chosen margin in sunlight or shadow, subject to reserve.
% Exported margin consumes battery energy; charging uses only solar above the target supply.
supplementalTarget_W = load_W * (1 + s.supplementalMargin_percent/100);
supplemental = comparisonBattery(max(supplementalTarget_W - solar_W, 0), ...
    max(solar_W - supplementalTarget_W, 0), dt, s);
result.comparison.labels = s.balanceLabels;
result.comparison.supply_W = [solar_W; batteryOnly.bus_W; ...
    solar_W + maximum.bus_W; solar_W + conditional.bus_W; solar_W + supplemental.bus_W];
result.comparison.balance_percent = 100 * (result.comparison.supply_W - load_W) ./ load_W;
result.comparison.batteryOnly = batteryOnly;
result.comparison.solarBatteryMaximum = maximum;
result.comparison.conditional = conditional;
result.comparison.supplemental = supplemental;
result.comparison.supplemental.targetMargin_percent = s.supplementalMargin_percent;
result.comparison.supplemental.targetSupply_W = supplementalTarget_W;
result.comparison.supplemental.targetShortfall_W = ...
    max(supplementalTarget_W - solar_W - supplemental.bus_W, 0);
result.comparison.note = ['Instantaneous supply margin; battery discharge uses stored energy. ' ...
    'Solar-Battery is a maximum-output test without recharge, not a load-following controller.'];

%% Battery discharge trade section
% Reuse the same solar/thermal solution and battery limits for independent discharge-limit cases.
% This is a one-orbit discharge experiment; exported surplus is not recycled into the battery.
if s.runBatteryTrade
    result.batteryDischargeTrade = finalPowerEPS.batteryDischargeTrade(t, solar_W, load_W, s);
    if s.showSummary
        trade = result.batteryDischargeTrade;
        if isempty(trade.bestIndex)
            fprintf('Battery trade: no sampled discharge limit meets +%g%% headroom and restores initial SOC.\n', ...
                trade.targetMargin_percent);
        else
            fprintf('Battery trade: lowest sampled limit for +%g%% headroom plus SOC recovery: %g%% (%.3f W pack).\n', ...
                trade.targetMargin_percent, trade.command_percent(trade.bestIndex), ...
                trade.commandPack_W(trade.bestIndex));
        end
    end
else
    result.batteryDischargeTrade = [];
end

%% Model checks section
% The operating range is a check, not a forced temperature clamp.
% Limits and charge efficiency remain placeholders until the array/battery ICDs are available.
if s.warnModelLimits && any(array.thermal.outsideProductOperatingRange)
    warning('powerBalance:productTemperature', ...
        'Settled temperatures [%g,%g] C leave the selected array range [%g,%g] C.', ...
        min(array.cellTemperature_C), max(array.cellTemperature_C), ...
        array.params.productOperatingRange_C(1), array.params.productOperatingRange_C(2));
end
if s.warnModelLimits && array.thermal.maxRadiationBiot > array.params.lumpedBiotWarningLimit
    warning('powerBalance:lumpedModel', 'The equivalent plate exceeds the radiation Biot screening limit.');
end
if s.showSummary
    fprintf('\nWorst-shadow envelope: %s\n', status);
    disp(result.summary)
    printTableKey(result.summary);
    fprintf('Worst-orbit average solar power (planning, includes eclipse): %.3f W.\n', ...
        result.solarPowerStats.orbitAveragePlanning_W);
    fprintf('Minimum solar power: %.3f W over the orbit; %.3f W during sunlit intervals.\n', ...
        result.solarPowerStats.minimumPlanning_W, result.solarPowerStats.minimumSunlitPlanning_W);
    fprintf('Required eclipse battery draw: %.3f Wh at the pack; peak delivered: %.2f W.\n', ...
        requiredPackDischarge_Wh, max(batteryBus_W));
    fprintf('Minimum charge efficiency without a rate limit: %.1f%%; assumed: %.1f%%.\n', ...
        100*minimumChargeEfficiency, 100*s.batteryChargeEfficiency);
    fprintf('Thermal orbit: %.2f to %.2f C; periodic mismatch %.4g C after %d repeats.\n', ...
        min(array.cellTemperature_C), max(array.cellTemperature_C), periodicError_C, orbit);
    fprintf('Published 35 W is per wing; this %.4f m^2 area model is not a rated-wing calibration.\n', ...
        s.arrayArea_m2);
end

%% Graphing section
% Five independent supply cases are normalized by load, with zero centered on the axis.
% Different line styles keep overlapping solar and combined histories distinguishable.
if s.makePlot
    figure('WindowStyle','normal','Color',s.figureColor,'Position',s.figurePosition,'Name','Worst-case orbit power balance');
    balance = result.comparison.balance_percent;
    limit = max(100, max(abs(balance), [], 'all')) * s.balanceAxisPadding;
    limit = ceil(limit / s.balanceTickStep_percent) * s.balanceTickStep_percent;
    patch([0 eclipse_s eclipse_s 0]/60, [-limit -limit limit limit], s.eclipseShadeColor, ...
        'EdgeColor','none','FaceAlpha',s.eclipseShadeAlpha,'HandleVisibility','off');
    hold on;
    yline(0, '-', 'Color',[0.30 0.30 0.30], 'HandleVisibility','off');
    curves = gobjects(5,1);
    for k = 1:5
        curves(k) = stairs(t/60, [balance(k,:) balance(k,end)], ...
            'Color',s.balanceColors(k,:), 'LineStyle',s.balanceLineStyles{k}, ...
            'LineWidth',s.lineWidth);
    end
    grid on; xlim([0 period_s]/60); ylim([-limit limit]);
    yticks(-limit:s.balanceTickStep_percent:limit); ytickformat('%g%%');
    xlabel('Time from eclipse entry [min]'); ylabel('Power balance relative to load [%]');
    title('Power balance');
    legend(curves, s.balanceLabels, 'Location','southoutside','NumColumns',2);
end

%% Battery trade graphing section
% Primary maps show power margin and SOC versus time and command.
% Optional four-view surfaces remain available; the black highlight is the lowest feasible sampled command.
if s.makePlot && s.runBatteryTrade && s.makeBatteryTradePlots
    plotDischargeTrade(result.batteryDischargeTrade, t, s);
end
end

%% Comparison battery function
% Track finite stored energy and conversion losses independently for each plotted case.
% Returned powers are interval averages; the last discharge interval may be energy-limited.
function history = comparisonBattery(requestedBus_W, chargeSurplus_W, dt, s)
capacity_Wh = s.batteryCapacityBOL_Wh * s.batteryCapacityRetained;
floor_Wh = capacity_Wh * s.batteryMinimumSOC;
ceiling_Wh = capacity_Wh * s.batteryMaximumSOC;
history.storedEnergy_Wh = zeros(1,numel(dt)+1);
history.storedEnergy_Wh(1) = capacity_Wh * s.batteryInitialSOC;
history.bus_W = zeros(size(dt));
history.packDischarge_W = zeros(size(dt));
history.packCharge_W = zeros(size(dt));
for k = 1:numel(dt)
    stored = history.storedEnergy_Wh(k);
    energyLimit_W = max(stored - floor_Wh,0) * 3600 / dt(k);
    discharge_W = min([requestedBus_W(k)/s.batteryToBusEfficiency, ...
        s.batteryMaxDischarge_W, energyLimit_W]);
    charge_W = 0;
    if requestedBus_W(k) == 0
        roomLimit_W = max(ceiling_Wh - stored,0) * 3600 / dt(k);
        charge_W = min([chargeSurplus_W(k)*s.batteryChargeEfficiency, ...
            s.batteryMaxCharge_W, roomLimit_W]);
    end
    history.packDischarge_W(k) = discharge_W;
    history.packCharge_W(k) = charge_W;
    history.bus_W(k) = discharge_W * s.batteryToBusEfficiency;
    history.storedEnergy_Wh(k+1) = stored + (charge_W-discharge_W)*dt(k)/3600;
end
history.SOC = history.storedEnergy_Wh / capacity_Wh;
end

%% Primary battery trade display
% Independent coordinates (time, commanded pack-power fraction) form a map.
% Flat faces preserve interval averages and abrupt transitions; no smooth curtain
% interpolates through an eclipse boundary or depletion event. Rate bins display
% the nearest sampled command, not a new interpolated operating point.
function plotDischargeTrade(trade,t,s)
rates=trade.command_percent(:); n=numel(rates);
if n==1
    rateEdges=[max(0,rates-0.5);min(100,rates+0.5)];
else
    rateEdges=[max(0,rates(1)-(rates(2)-rates(1))/2); ...
        (rates(1:end-1)+rates(2:end))/2;min(100,rates(end)+(rates(end)-rates(end-1))/2)];
end
fig=figure('WindowStyle','normal','Color',s.figureColor, ...
    'Name','Battery limit trade');
finalPowerEPS.fitFigureToScreen(fig,s.batteryTradeFigureSize_px);
layout=tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');
a=nexttile(layout); mapCells(a,t/60,rateEdges,trade.balance_percent);
limit=max(100,max(abs(trade.balance_percent(:))));
colormap(a,divergingMap(257)); caxis(a,[-limit limit]);
cb=colorbar(a); cb.Label.String='Available margin [% of load]';
title(a,'Available power margin');
limitLabel=sprintf('Battery power limit [%% of %g W]',s.batteryMaxDischarge_W);
ylabel(a,limitLabel); hold(a,'on');
mid=(t(1:end-1)+t(2:end))/120;
mask=double(trade.balance_percent>=trade.targetMargin_percent-s.batteryTradeTolerance_percent);
if n>1 && numel(mid)>1 && any(mask(:)==0) && any(mask(:)==1)
    contour(a,mid,rates,mask,[0.5 0.5],'k:','LineWidth',1.2,'HandleVisibility','off');
end
b=nexttile(layout); mapCells(b,t/60,rateEdges,100*trade.SOC(:,2:end));
colormap(b,parula(256)); caxis(b,100*[s.batteryMinimumSOC s.batteryMaximumSOC]);
cb=colorbar(b); cb.Label.String='SOC at interval end [%]';
title(b,'Battery state of charge'); xlabel(b,'Time from eclipse entry [min]');
ylabel(b,limitLabel);
for ax=[a b]
    hold(ax,'on'); xline(ax,s.eclipseDuration_min,'k--','HandleVisibility','off');
    if ~isempty(trade.bestIndex)
        yline(ax,rates(trade.bestIndex),'k-','LineWidth',2.2,'HandleVisibility','off');
    end
end
h=plot(a,NaN,NaN,'k--'); names={'Eclipse exit'};
h(end+1)=plot(a,NaN,NaN,'k:'); names{end+1}=sprintf('%g%% margin boundary',trade.targetMargin_percent);
if ~isempty(trade.bestIndex)
    h(end+1)=plot(a,NaN,NaN,'k-','LineWidth',2.2); names{end+1}=sprintf('Lowest feasible limit: %g%%',rates(trade.bestIndex));
end
lg=legend(a,h,names,'Orientation','horizontal','NumColumns',3,'FontSize',9,'Box','off'); lg.Layout.Tile='south';
set([a b],'FontSize',10);
linkaxes([a b],'x');

% Three synchronized histories distinguish capability, actual output, and charge.
targets=s.batteryTradeCurveCommands_percent(:); rows=zeros(size(targets));
for k=1:numel(targets), [~,rows(k)]=min(abs(rates-targets(k))); end
rows=unique([rows;trade.bestIndex]); colors=lines(numel(rows));
fig=figure('WindowStyle','normal','Color',s.figureColor, ...
    'Name','Battery operation');
finalPowerEPS.fitFigureToScreen(fig,s.batteryOperationFigureSize_px);
layout=tiledlayout(fig,3,1,'TileSpacing','compact','Padding','compact');
a=nexttile(layout); hold(a,'on'); b=nexttile(layout); hold(b,'on'); c=nexttile(layout); hold(c,'on');
h=gobjects(numel(rows),1); names=cell(numel(rows),1);
% Draw the highlighted case last so coincident higher-limit traces cannot hide it.
if isempty(trade.bestIndex)
    plotOrder=(1:numel(rows)).';
else
    plotOrder=[find(rows~=trade.bestIndex);find(rows==trade.bestIndex)];
end
styles={'-','--',':','-.'};
for j=plotOrder(:).'
    k=rows(j); color=colors(j,:); width=1.8; style=styles{1+mod(j-1,numel(styles))};
    names{j}=sprintf('%g%% limit (%.3g W)',rates(k),trade.commandPack_W(k));
    if isequal(k,trade.bestIndex), color=[0 0 0]; width=2.6; style='-'; names{j}=[names{j} ' (min.)']; end
    h(j)=stairs(a,t/60,[trade.balance_percent(k,:) trade.balance_percent(k,end)],'Color',color,'LineWidth',width,'LineStyle',style);
    stairs(b,t/60,[trade.packDischarge_W(k,:) trade.packDischarge_W(k,end)],'Color',color,'LineWidth',width,'LineStyle',style);
    plot(c,t/60,100*trade.SOC(k,:),'Color',color,'LineWidth',width,'LineStyle',style);
end
% Labels belong in the external legend; no words cover the curves.
yline(a,0,'k-','HandleVisibility','off'); yline(a,trade.targetMargin_percent,'k:','HandleVisibility','off');
yline(c,100*s.batteryMinimumSOC,'k:','HandleVisibility','off');
for ax=[a b c]
    grid(ax,'on'); xlim(ax,[t(1) t(end)]/60);
    xline(ax,s.eclipseDuration_min,'k--','HandleVisibility','off');
end
h(end+1)=plot(a,NaN,NaN,'k--'); names{end+1}='Eclipse exit';
h(end+1)=plot(a,NaN,NaN,'k:'); names{end+1}=sprintf('Target %g%% / reserve %g%%',trade.targetMargin_percent,100*s.batteryMinimumSOC);
title(a,'Available power margin'); ylabel(a,'Available margin [%]'); ytickformat(a,'%g%%');
title(b,'Actual discharge'); ylabel(b,'Pack power [W]');
title(c,'Battery state of charge'); ylabel(c,'SOC [%]'); xlabel(c,'Time from eclipse entry [min]');
upper=max(trade.packDischarge_W(rows,:),[],'all');
ylim(b,[0 max(1,1.12*upper)]);
margin=trade.balance_percent(rows,:); low=min(0,min(margin(:))); high=max(0,max(margin(:)));
padding=0.05*max(1,high-low); ylim(a,[low-padding high+padding]);
ylim(c,[100*s.batteryMinimumSOC-2 100*s.batteryMaximumSOC+2]);
lg=legend(a,h,names,'Orientation','horizontal','NumColumns',3,'FontSize',9,'Box','off'); lg.Layout.Tile='south';
set([a b c],'FontSize',10);
linkaxes([a b c],'x');
if s.makeBatteryTradeSurfaceViews, plotDischargeTradeSurface(trade,t,s); end
end
function mapCells(ax,xEdges,yEdges,values)
% pcolor omits its last vertex row/column: append those vertices deliberately.
C=[values values(:,end)]; C=[C;C(end,:)];
pcolor(ax,xEdges,yEdges,C); shading(ax,'flat'); set(ax,'YDir','normal');
xlim(ax,[xEdges(1) xEdges(end)]); ylim(ax,[yEdges(1) yEdges(end)]);
end
function colors=divergingMap(n)
colors=interp1([0 0.5 1],[0.75 0.10 0.12;1 1 1;0.05 0.35 0.75],linspace(0,1,n));
end

%% Battery trade graphing function
% X=time, Y=power balance relative to load, Z=commanded percent of pack limit, color=SOC.
% Balance = 100*(solar bus power + available battery capability - load)/load.
% Zero meets the load, positive is available headroom, negative is a shortfall. The battery
% command stays on Z even when the reserve floor limits actual power. Balance uses
% interval averages; SOC is sampled at interval boundaries, as in the original surface.
function plotDischargeTradeSurface(trade, t, s)
colors = interp1(linspace(0,1,size(s.batteryTradeColorAnchors,1)), ...
    s.batteryTradeColorAnchors, linspace(0,1,s.batteryTradeColorLevels));
[timeGrid,commandGrid] = meshgrid(t/60,trade.command_percent);
balance = [trade.balance_percent trade.balance_percent(:,end)];
socPercent = 100*trade.SOC;
balanceLimit = max(100,max(abs(balance),[],'all')) * s.balanceAxisPadding;
balanceLimit = ceil(balanceLimit/s.balanceTickStep_percent) * s.balanceTickStep_percent;
views = [s.batteryTradeView_deg(:).'; 0 90; 0 0; 90 0];
viewNames = {'3D','XY','XZ','YZ'};

% Select the same sparse mesh once, then reuse it in each camera view.
rateRows = []; timeCols = [];
if numel(trade.command_percent) > 1
    rateTargets = unique([trade.command_percent(1):s.batteryTradeMeshRateStep_percent: ...
        trade.command_percent(end), trade.command_percent(end)]);
    rateRows = zeros(size(rateTargets));
    for j = 1:numel(rateTargets)
        [~,rateRows(j)] = min(abs(trade.command_percent-rateTargets(j)));
    end
    rateRows = unique(rateRows);
    timeTargets = unique([t(1)/60:s.batteryTradeMeshTimeStep_min:t(end)/60, t(end)/60]);
    timeCols = zeros(size(timeTargets));
    for j = 1:numel(timeTargets)
        [~,timeCols(j)] = min(abs(t/60-timeTargets(j)));
    end
    timeCols = unique(timeCols);
end
for v = 1:numel(viewNames)
    fig = figure('WindowStyle','normal','Color',s.figureColor,'Position',s.figurePosition, ...
        'Name',['Battery discharge trade: ' viewNames{v} ' view']);
    ax = axes('Parent',fig);
    if numel(trade.command_percent) > 1
        surf(ax,timeGrid,balance,commandGrid,socPercent, ...
            'EdgeColor','none','FaceColor','interp','FaceAlpha',s.batteryTradeSurfaceAlpha);
    else
        % A duplicated row gives a color-interpolated 3D line for a one-rate trade.
        surface(ax,[timeGrid;timeGrid],[balance;balance],[commandGrid;commandGrid], ...
            [socPercent;socPercent],'FaceColor','none','EdgeColor','interp', ...
            'LineWidth',s.batteryTradeLineWidth);
    end
    hold(ax,'on');
    for k = rateRows
        plot3(ax,timeGrid(k,:),balance(k,:),commandGrid(k,:), ...
            'Color',s.batteryTradeMeshColor,'LineWidth',s.batteryTradeMeshLineWidth, ...
            'HandleVisibility','off');
    end
    for j = timeCols
        plot3(ax,timeGrid(:,j),balance(:,j),commandGrid(:,j), ...
            'Color',s.batteryTradeMeshColor,'LineWidth',s.batteryTradeMeshLineWidth, ...
            'HandleVisibility','off');
    end
    if ~isempty(trade.bestIndex)
        k = trade.bestIndex;
        best = plot3(ax,timeGrid(k,:),balance(k,:),commandGrid(k,:), 'k-', ...
            'LineWidth',s.batteryTradeHighlightWidth);
        legend(ax,best,sprintf('Lowest feasible sample: %g%%',trade.command_percent(k)), ...
            'Location','southoutside');
    end
    grid(ax,'on'); view(ax,views(v,:)); xlim(ax,[t(1) t(end)]/60);
    ylim(ax,[-balanceLimit balanceLimit]); zlim(ax,[0 100]);
    yticks(ax,-balanceLimit:s.balanceTickStep_percent:balanceLimit); ytickformat(ax,'%g%%');
    xlabel(ax,'Time from eclipse entry [min]'); ylabel(ax,'Available headroom [% of load]');
    zlabel(ax,'Discharge limit [% of pack maximum]');
    title(ax,['Battery trade: ' viewNames{v}]);
    colormap(ax,colors); caxis(ax,[0 100]);
    cb = colorbar(ax); cb.Label.String = 'Battery state of charge [%]';
end
end



%% Table key function
% Print definitions stored with each table so column names and explanations stay together.
function printTableKey(values)
fprintf('Table key (BOL = beginning of life; planning = after aging allowance):\n');
for k = 1:width(values)
    fprintf('  %-29s %s\n', values.Properties.VariableNames{k}, ...
        values.Properties.VariableDescriptions{k});
end
fprintf('\n');
end

%% References
% [NA] Kongsberg NanoAvionics, CubeSat GaAs Solar Panels; checked 2026-10-19.
% https://nanoavionics.com/cubesat-components/cubesat-gaas-solar-panel/
% Product inputs, thermal proxies and their sources are listed in arrayPower.m.
%
% [BT] Team battery architecture trade, uploaded "Battery_Trade 3.pdf", Rev. 1.
% Section 1.02: 550 km / 35.8 min shadow, eta_dis=0.95, SOH=0.80, 40% reserve floor.
% Section 1.07: Candidate A 45 Wh and 50 W peak-power estimate; 1.08 selects that architecture.
% The new 23.65 W continuous load comes from the user's instruction, replacing the mode schedule.
% A final battery model, continuous rating, voltage dependence and charging limits are not selected.
%
% [ENV] NASA/TM-2001-211222, Simple Thermal Environment Model, Sections 1.2.2-1.2.3.
% https://ntrs.nasa.gov/api/citations/20020001443/downloads/20020001443.pdf
% Uses Earth orbital eccentricity 0.0167 for solar-distance extrema; here S=cfg.S0/1.0167^2.
% The model retains the project's cfg.S0, not the older report's nominal solar constant.
%
% [ENERGY] Conservation of energy: dE/dt=eta_charge*P_surplus-P_battery_bus/eta_discharge.
% Orbital period follows T=2*pi*sqrt(a^3/mu), as in the team's propagate.m / Curtis reference.
% arrayPower now uses front/rear thermal control volumes; incidence cases are independent orientations, not spatial nodes.
