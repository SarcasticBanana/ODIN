function h=solarArrayElectricalPower(irradiance_W_m2,temperature_K,p)
% SOLARARRAYELECTRICALPOWER Shared original temperature/power equation.
% Retains all existing efficiencies, coverage, loss and aging assumptions.
active=irradiance_W_m2>0 & p.selfShadowRetained>0 & p.cellCoverage>0 & p.cellEfficiency>0;
if isscalar(irradiance_W_m2) && isscalar(temperature_K)
    factor=1;
    if active, factor=1+p.powerTempCoeff_perC*(temperature_K-273.15-p.referenceTemperature_C); end
    h.cellEfficiency=p.cellEfficiency*factor;
else
    factor=ones(size(irradiance_W_m2+temperature_K));
    correction=1+p.powerTempCoeff_perC.*(temperature_K-273.15-p.referenceTemperature_C);
    correction=correction+zeros(size(factor)); active=active & true(size(factor));
    factor(active)=correction(active);
    h.cellEfficiency=p.cellEfficiency.*factor;
end
if any(h.cellEfficiency(active)<0 | h.cellEfficiency(active)>1)
    error('solarArrayHeatBalance:efficiency','Temperature extrapolation implies an invalid active cell efficiency.');
end
retained=1;
if isfield(p,'electricalDegradationRetained'), retained=p.electricalDegradationRetained; end
h.absorbedSolar_W_m2=p.solarAbsorptivity.*irradiance_W_m2.*p.selfShadowRetained;
h.potentialElectrical_W_m2=irradiance_W_m2.*h.cellEfficiency.*p.cellCoverage.* ...
    p.selfShadowRetained.*p.mismatchRetained.*p.opticalRetained.*retained;
if any(h.potentialElectrical_W_m2(:)>h.absorbedSolar_W_m2(:)+p.energyBalanceTolerance_W_m2)
    error('solarArrayHeatBalance:energy','Available electrical power exceeds absorbed sunlight.');
end
end
