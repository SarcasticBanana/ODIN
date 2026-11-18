function testComponentBoundaries()
% TESTCOMPONENTBOUNDARIES Independent closed-form boundary and sign checks.
q=struct('showSummary',false,'warnAssumptions',false);
b=finalPowerEPS.BatteryThermalModel(struct('packCurrent_A',[2 -2],'cellVoltage_V',3.7,'time_s',[0 1]),q);
assert(max(abs(b.Qirreversible_W-.36))<1e-12); % 2S2P, 1 A/cell, 0.09 ohm/cell
assert(max(abs(b.packTerminalPower_W-[14.8 -14.8]))<1e-12);
in=struct('packCurrent_A',2,'cellVoltage_V',3.6,'cellOCVoltage_V',3.7);
b=finalPowerEPS.BatteryThermalModel(in,q); assert(abs(b.Qirreversible_W-.4)<1e-12); % no additional I²R
q2=q; q2.entropyCoefficient_V_K=.001;
b=finalPowerEPS.BatteryThermalModel(struct('packCurrent_A',2,'cellTemperature_C',26.85),q2);
assert(abs(b.Qreversible_W+1.2)<1e-12);
q2=q; q2.packInterconnectResistance_Ohm=.05;
b=finalPowerEPS.BatteryThermalModel(struct('packTerminalPower_W',10),q2);
assert(abs(b.packTerminalPower_W-10)<1e-12);
m=finalPowerEPS.PMADThermalModel(struct('solarToLoad_W',10,'batteryToLoad_W',0,'batteryCharge_W',5),q);
assert(abs(m.Qgen_W-15*(1/.931-1))<1e-12);
assert(abs(m.totalElectricalInput_W-m.totalElectricalOutput_W-m.Qgen_W)<1e-12);
rejected=false;
try
 finalPowerEPS.PMADThermalModel(struct('solarToLoad_W',0,'batteryToLoad_W',1,'batteryCharge_W',1),q);
catch err
 rejected=strcmp(err.identifier,'PMADThermalModel:batteryFlow');
end
assert(rejected);
rejected=false;
try
 finalPowerEPS.SolarArrayThermalModel(struct('cellTemperature_C',20,'temperatureRate_K_s',0),q);
catch err
 rejected=strcmp(err.identifier,'SolarArrayThermalModel:temperature');
end
assert(rejected);
a=finalPowerEPS.SolarArrayThermalModel(struct('irradiance_W_m2',1000,'cellTemperature_C',28,'electricalPower_W',100),q);
assert(a.extractionExceedsAvailable); % absorption bound alone does not prove electrical feasibility
fprintf('PASS: 10 component boundary, sign, infeasibility and derivative checks.\n');
end
