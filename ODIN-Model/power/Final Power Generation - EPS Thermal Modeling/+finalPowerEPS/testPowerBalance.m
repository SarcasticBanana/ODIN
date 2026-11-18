function testPowerBalance()
% TESTPOWERBALANCE Focused checks for energy conservation and unavailable-power cases.
% Run in MATLAB with the updated files and the team's config.m / solarFlux.m on the path.

%% Default settings checks
p = finalPowerEPS.arrayPower();
assert(p.cellEfficiency == 0.295);
assert(abs(p.panelThickness_m - 0.0013) < 1e-12);
assert(p.selectedArrayArea_m2 > 0);
assert(p.backEmissivity > 0 && p.thermalNodeCount == 2);
cfg = finalPowerEPS.config();
settings = struct('makePlot',false,'showSummary',false,'warnModelLimits',false);

%% Nominal balance checks
r = finalPowerEPS.powerBalance(cfg, settings);
assert(r.params.arrayArea_m2 == p.selectedArrayArea_m2);
dt_h = r.intervalDuration_s / 3600;
assert(r.loadCovered && r.energyPositive);
assert(all(r.batteryToBus_W(r.sunlit) == 0));
assert(all(r.solarAvailable_W(~r.sunlit) == 0));
assert(abs(r.energy.load_Wh - 23.65*r.orbitPeriod_min/60) < 1e-9);
assert(abs(r.energy.requiredPackDischarge_Wh - 23.65*35.8/60/0.95) < 1e-9);
assert(max(abs(diff(r.storedEnergy_Wh) - ...
    (r.batteryPackCharge_W-r.batteryPackDischarge_W).*dt_h)) < 1e-10);
assert(max(abs(min(r.solarAvailable_W,r.load_W) + ...
    r.batteryToBus_W + r.unserved_W - r.load_W)) < 1e-9);
assert(min(r.batterySOC) >= r.params.batteryMinimumSOC - 1e-10);
assert(max(r.batterySOC) <= r.params.batteryMaximumSOC + 1e-10);
assert(abs(r.batterySOC(end)-r.batterySOC(1)) < 1e-8);

%% Insufficient power and storage checks
settings.batteryMaxDischarge_W = 10;
weak = finalPowerEPS.powerBalance(cfg, settings);
assert(~weak.loadCovered && weak.energy.unserved_Wh > 0);
assert(max(weak.batteryToBus_W) <= 9.5 + 1e-9);
settings = rmfield(settings, 'batteryMaxDischarge_W');
settings.batteryCapacityBOL_Wh = 10;
small = finalPowerEPS.powerBalance(cfg, settings);
assert(~small.loadCovered && small.energy.unserved_Wh > 0);
assert(min(small.batterySOC) >= small.params.batteryMinimumSOC - 1e-10);
settings = rmfield(settings, 'batteryCapacityBOL_Wh');

%% Edge-on projection and charging checks
settings.arrayOverrides = struct('areaStudyIncidence_deg',90);
edge = finalPowerEPS.powerBalance(cfg, settings);
assert(all(edge.solarAvailable_W == 0)); % cosd(90), without an angle-specific output switch
assert(~edge.energyPositive && ~edge.loadCovered);
settings = rmfield(settings, 'arrayOverrides');
settings.batteryChargeEfficiency = 0;
blocked = finalPowerEPS.powerBalance(cfg, settings);
assert(~blocked.energyPositive && all(blocked.batteryPackCharge_W == 0));
settings = rmfield(settings, 'batteryChargeEfficiency');
settings.batteryMaxCharge_W = 1;
slow = finalPowerEPS.powerBalance(cfg, settings);
assert(~slow.energyPositive);
fprintf('Power-balance conservation and limiting-case checks passed.\n');
end
