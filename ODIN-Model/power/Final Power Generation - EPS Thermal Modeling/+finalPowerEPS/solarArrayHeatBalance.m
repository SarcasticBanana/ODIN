function h = solarArrayHeatBalance(irradiance_W_m2, temperature_K, p, electrical_W_m2, heater_W_m2, back_K, diagnostics)
% SOLARARRAYHEATBALANCE Shared panel-boundary electrical/thermal equations.
% Internal kernel for arrayPower and SolarArrayThermalModel; inputs are validated by callers.
% All heat/power outputs are W per m^2 of physical front-face area; positive heat flows OUT.
% Empty electrical_W_m2 requests the available array-terminal power at this temperature.
% Controller/MPPT/harness conversion losses are outside this panel boundary.
if nargin < 4, electrical_W_m2 = []; end
if nargin < 5, heater_W_m2 = 0; end
if nargin < 6 || isempty(back_K), back_K=temperature_K; end
if nargin<7, diagnostics=true; end
nodes=1; if isfield(p,'thermalNodeCount'), nodes=p.thermalNodeCount; end
sigma = 5.670374419e-8;
h=finalPowerEPS.solarArrayElectricalPower(irradiance_W_m2,temperature_K,p);
h.frontAbsorbedEarth_W_m2=p.absorbedEarthFlux_W_m2+p.frontAbsorbedEarthFlux_W_m2;
h.backAbsorbedEarth_W_m2=p.backAbsorbedEarthFlux_W_m2;
h.absorbedEarth_W_m2=h.frontAbsorbedEarth_W_m2+h.backAbsorbedEarth_W_m2+zeros(size(temperature_K));
if isempty(electrical_W_m2), electrical_W_m2=h.potentialElectrical_W_m2; end
h.electrical_W_m2 = electrical_W_m2 + zeros(size(irradiance_W_m2+temperature_K));
h.heater_W_m2 = heater_W_m2 + zeros(size(irradiance_W_m2+temperature_K));
h.heatInput_W_m2 = h.absorbedSolar_W_m2 + h.absorbedEarth_W_m2 - h.electrical_W_m2 + h.heater_W_m2;
frontSpace=sigma*p.frontEmissivity.*(1-p.frontStructureViewFactor-p.frontEarthViewFactor).*(temperature_K.^4-p.spaceTemperature_K^4);
backSpace=sigma*p.backEmissivity.*(1-p.backStructureViewFactor-p.backEarthViewFactor).*(back_K.^4-p.spaceTemperature_K^4);
frontStructure=sigma*p.frontEmissivity*p.frontStructureViewFactor.*(temperature_K.^4-(p.structureTemperature_C+273.15).^4);
backStructure=sigma*p.backEmissivity*p.backStructureViewFactor.*(back_K.^4-(p.structureTemperature_C+273.15).^4);
frontEarth=sigma*p.frontEmissivity.*p.frontEarthViewFactor.*temperature_K.^4;
backEarth=sigma*p.backEmissivity.*p.backEarthViewFactor.*back_K.^4;
h.radiationToEarth_W_m2=frontEarth+backEarth;
h.frontRadiation_W_m2=frontSpace+frontStructure+frontEarth;
h.backRadiation_W_m2=backSpace+backStructure+backEarth;
h.radiationToSpace_W_m2=frontSpace+backSpace;
h.radiationToStructure_W_m2=frontStructure+backStructure;
h.conductionToStructure_W_m2 = p.mountConductancePerArea_W_m2K .* ...
    ((temperature_K+back_K)/2-273.15-p.mountTemperature_C);
h.heatToStructure_W_m2 = h.radiationToStructure_W_m2 + h.conductionToStructure_W_m2;
h.heatOut_W_m2 = h.radiationToSpace_W_m2 + h.radiationToEarth_W_m2 + h.heatToStructure_W_m2;
if diagnostics
coefficient=@(value,sink) sigma.*(value+sink).*(value.^2+sink.^2);
h.radiationCoefficient_W_m2K=p.frontEmissivity.*((1-p.frontStructureViewFactor-p.frontEarthViewFactor).*coefficient(temperature_K,p.spaceTemperature_K)+ ...
    p.frontStructureViewFactor.*coefficient(temperature_K,p.structureTemperature_C+273.15)+p.frontEarthViewFactor.*coefficient(temperature_K,0)) + ...
    p.backEmissivity.*((1-p.backStructureViewFactor-p.backEarthViewFactor).*coefficient(back_K,p.spaceTemperature_K)+ ...
    p.backStructureViewFactor.*coefficient(back_K,p.structureTemperature_C+273.15)+p.backEarthViewFactor.*coefficient(back_K,0));
end
h.storage_W_m2 = h.heatInput_W_m2-h.heatOut_W_m2;
h.arealHeatCapacity_J_m2K = p.panelDensity_kg_m3*p.panelThickness_m*p.panelSpecificHeat_J_kgK + ...
    p.additionalArealHeatCapacity_J_m2K;
h.frontToBack_W_m2=zeros(size(h.storage_W_m2));
if nodes==2
    % Half-thickness control volumes; their centroids are d/2 apart. G/A=2*k/d.
    h.frontToBack_W_m2=2*p.panelConductivity_W_mK/p.panelThickness_m.*(temperature_K-back_K);
    frontMount=0.5*p.mountConductancePerArea_W_m2K.*(temperature_K-273.15-p.mountTemperature_C);
    backMount=0.5*p.mountConductancePerArea_W_m2K.*(back_K-273.15-p.mountTemperature_C);
    h.temperatureRate_K_s=2*(h.heatInput_W_m2-h.backAbsorbedEarth_W_m2-h.frontRadiation_W_m2-frontMount-h.frontToBack_W_m2)/h.arealHeatCapacity_J_m2K;
    h.backTemperatureRate_K_s=2*(h.backAbsorbedEarth_W_m2+h.frontToBack_W_m2-h.backRadiation_W_m2-backMount)/h.arealHeatCapacity_J_m2K;
else
    h.temperatureRate_K_s=h.storage_W_m2/h.arealHeatCapacity_J_m2K;
    h.backTemperatureRate_K_s=h.temperatureRate_K_s;
end
end
