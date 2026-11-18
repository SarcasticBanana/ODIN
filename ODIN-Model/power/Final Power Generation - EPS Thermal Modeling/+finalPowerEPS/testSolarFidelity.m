function testSolarFidelity(runOrbit)
% TESTSOLARFIDELITY Analytic physics, I-V calibration, and integrated orbit checks.
% Run with this package and the team's orbital dependencies on the MATLAB path.
if nargin==0, runOrbit=true; end
count=0; quiet=struct('showSummary',false,'warnAssumptions',false);
q=quiet; q.area_m2=.175; q.initialTemperature_C=20; q.initialBackTemperature_C=20;
q.frontEmissivity=0; q.backEmissivity=0; q.defaultHeaterPower_W=10;
x=finalPowerEPS.SolarArrayThermalModel(struct('timeEdges_s',0:10:600,'irradiance_W_m2',0),q);
C=x.thermalCapacitance_J_K; p=x.params; g=2*p.panelConductivity_W_mK/p.panelThickness_m;
CA=C/p.area_m2; flux=10/p.area_m2;
meanExact=20+10*x.timeEdges_s/C;
deltaExact=flux/(2*g)*(1-exp(-4*g*x.timeEdges_s/CA));
near((x.temperatureEdges_C+x.backTemperatureEdges_C)/2,meanExact,1e-8,'Total mass/heat capacity is counted once.');
near(x.temperatureEdges_C-x.backTemperatureEdges_C,deltaExact,2e-5,'Two-node conduction matches the analytic differential mode.');
near(x.QheatInput_W,10+zeros(size(x.time_s)),1e-10,'Heater input is counted once.');
near(x.balanceResidual_W,zeros(size(x.time_s)),1e-7,'Node exchange cancels in the total energy balance.');
q=quiet; q.spaceTemperature_K=0; q.initialTemperature_C=80; q.initialBackTemperature_C=80;
x=finalPowerEPS.SolarArrayThermalModel(struct('timeEdges_s',0:10:600,'irradiance_W_m2',0),q);
CA=x.thermalCapacitance_J_K/x.params.area_m2;
exact=((80+273.15)^-3+3*(x.params.frontEmissivity+x.params.backEmissivity)*5.670374419e-8*x.timeEdges_s/CA).^(-1/3)-273.15;
near(x.temperatureEdges_C,exact,2e-4,'Equal-temperature two-face vacuum cooling matches Stefan-Boltzmann analytic cooling.');
near(x.backTemperatureEdges_C,x.temperatureEdges_C,1e-8,'Symmetric dark faces stay at equal temperatures.');
near(sum(x.Qout_W.*diff(x.timeEdges_s)),x.thermalCapacitance_J_K*(80-x.temperatureEdges_C(end)),1e-6,'Radiated energy equals lost stored energy.');
q=quiet; q.frontStructureViewFactor=.2; q.backStructureViewFactor=.4; q.mountConductancePerArea_W_m2K=3;
x=finalPowerEPS.SolarArrayThermalModel(struct('cellTemperature_C',80,'backTemperature_C',60,'irradiance_W_m2',0),q);
near(x.QconductionToStructure_W,x.params.area_m2*3*(70-x.params.mountTemperature_C),1e-10,'Mount conductance is shared between nodes and counted once.');
near(x.Qout_W,x.QradiationToSpace_W+x.QtoStructure_W,1e-10,'All enabled outward paths sum to Qout.');
near(x.Qradiation_W,x.QfrontRadiation_W+x.QbackRadiation_W,1e-10,'Face-resolved radiation sums to panel radiation.');
q.structureTemperature_C=150; q.mountTemperature_C=150;
x=finalPowerEPS.SolarArrayThermalModel(struct('cellTemperature_C',80,'backTemperature_C',60,'irradiance_W_m2',0),q);
check(x.QtoStructure_W<0,'A hotter structure heats the panel with the correct sign.');
q=quiet; q.thermalNodeCount=1; q.backEmissivity=0; q.spaceTemperature_K=0;
x=finalPowerEPS.SolarArrayThermalModel(struct('timeEdges_s',0:10:600,'irradiance_W_m2',0),q);
CA=x.thermalCapacitance_J_K/x.params.area_m2;
exact=((x.params.initialTemperature_C+273.15)^-3+3*x.params.frontEmissivity*5.670374419e-8*x.timeEdges_s/CA).^(-1/3)-273.15;
near(x.temperatureEdges_C,exact,2e-4,'Legacy one-node/one-face option retains its physical cooling behavior.');

ivQuiet=struct('makePlots',false,'showSummary',false,'mode','xtj_reference');
g=finalPowerEPS.arrayGeometry(); area=g.cellArea_cm2;
in=struct('time_s',0,'irradiance_W_m2',1353,'cellTemperature_C',28);
r=finalPowerEPS.solarArrayIV(in,ivQuiet);
near(r.Voc_V,16*2.720,1e-8,'Series-cell count determines proxy Voc.');
near(r.Isc_A,4*.018*area,1e-9,'Current density times cell area and parallel strings gives Isc.');
near(r.Vmp_V,16*2.406,1e-7,'MPP voltage follows assumed integer series count.');
near(r.Imp_A,4*.0175*area,1e-9,'MPP current follows active area and parallel strings.');
check(r.referenceFit.maxRelativeResidual<1e-8,'Diode fit meets source endpoints and MPP slope.');
check(all(diff(r.curveCurrent_A)<=1e-12),'Forward I-V curve is monotone.');
check(max(r.curveCurrent_A.*r.curveVoltage_V)<=r.Pmp_W+1e-8,'MPP is the global forward power maximum.');
settings=ivQuiet; settings.testVoltage_V=r.Vmp_V-1e-4; left=finalPowerEPS.solarArrayIV(in,settings);
settings.testVoltage_V=r.Vmp_V+1e-4; right=finalPowerEPS.solarArrayIV(in,settings);
near((right.powerAtTestVoltage_W-left.powerAtTestVoltage_W)/2e-4,0,1e-6,'Independent finite difference verifies MPP.');
check(isempty(r.testVoltage_V),'No arbitrary fixed voltage is imposed.');
settings=ivQuiet; settings.testVoltage_V=0; shorted=finalPowerEPS.solarArrayIV(in,settings);
near(shorted.currentAtTestVoltage_A,r.Isc_A,1e-9,'Short-circuit test matches Isc.');
settings.testVoltage_V=r.Voc_V+1; opened=finalPowerEPS.solarArrayIV(in,settings);
check(opened.currentAtTestVoltage_A==0,'No forward generation above Voc.');
in.time_s=[0 1 2]; in.irradiance_W_m2=[1353 676.5 0];
r=finalPowerEPS.solarArrayIV(in,ivQuiet);
near(r.Isc_A(2)/r.Isc_A(1),.5,1e-6,'Photocurrent scales with irradiance.');
check(r.Voc_V(2)<r.Voc_V(1) && r.Voc_V(2)>.9*r.Voc_V(1),'Voc changes nonlinearly with illumination.');
check(r.Pmp_W(3)==0,'Eclipse has zero available forward generation.');
in=struct('time_s',[0 1 2],'irradiance_W_m2',1353,'cellTemperature_C',[15 28 85]);
r=finalPowerEPS.solarArrayIV(in,ivQuiet);
check(all(diff(r.Voc_V)<0)&all(diff(r.Vmp_V)<0),'Warming lowers voltage.');
near(r.Voc_V,16*(2.720-.0056*([15 28 85]-28)),1e-8,'Published Voc temperature slope is reproduced.');
check(~r.verifiedForNanoAvionics && ~r.wiringVerified,'No unsupported vendor/wiring verification claim.');
settings=ivQuiet; settings.seriesCells=8; eight=finalPowerEPS.solarArrayIV(in,settings);
near(eight.Vmp_V,r.Vmp_V/2,1e-8,'Half series count halves voltage for fixed 64-cell count.');
near(eight.Imp_A,2*r.Imp_A,1e-8,'More parallel strings increase current.');
near(eight.Pmp_W,r.Pmp_W,1e-8,'Rewiring fixed total cells preserves ideal total MPP power.');
in.availablePower_W=[40 30 0];
r=finalPowerEPS.solarArrayIV(in,struct('makePlots',false,'showSummary',false));
near(r.Pmp_W,in.availablePower_W,1e-8,'Proxy preserves shared panel power at every sample.');
check(all(r.currentNormalization>=0),'Current normalization is explicit and nonnegative.');
check(g.oldCellsTotalArea_m2>.175 && abs(64*g.cellArea_m2-g.totalCellArea_m2)<1e-12,'Old cells fail area constraint; resized footprints satisfy it.');
near(g.cellCoverage,64*.00310975/.1999,1e-12,'Coverage derives from the old CAD inputs.');
near(g.uniformLinearScale,sqrt(.175/.1999),1e-12,'Equal dimension scaling follows square root of area ratio.');

if runOrbit
    cfg=finalPowerEPS.config(); cfg.h=600; out=finalPowerEPS.propagate(cfg,1);
    settings=struct('runIncidenceSweep',true,'incidenceAngles_deg',[0 45 90], ...
        'incidenceSummaryAngles_deg',[0 45 90],'makeMainArrayPlots',false,'showAreaSummary',false, ...
        'warnTemperatureExtrapolation',false,'warnLumpedModel',false,'warnProductTemperature',false);
    a=finalPowerEPS.arrayPower(out,cfg,[],settings); x=finalPowerEPS.SolarArrayThermalModel(a,quiet);
    near(x.Qout_W,a.thermal.heatOut_W_m2*a.params.selectedArrayArea_m2,1e-8,'Array sizing and thermal postprocessing use identical face-resolved physics.');
    near(x.balanceResidual_W,zeros(size(x.time_s)),1e-7,'Independent front/rear ODE derivatives close the postprocessed balance.');
    check(all(a.incidenceSweep.arrayBOL_W(3,:)==0),'Edge-on illumination produces zero electricity.');
    settings.runIncidenceSweep=false; settings.backEmissivity=0;
    one=finalPowerEPS.arrayPower(out,cfg,[],settings);
    check(all(a.cellTemperature_C<=one.cellTemperature_C+1e-3),'Exposed rear radiation cools the front relative to the same two-node model without rear emission.');
    check(trapz(a.time_s,a.busPlanning_W(1,:))>trapz(one.time_s,one.busPlanning_W(1,:)),'Cooler cells increase generated energy under the retained negative coefficient.');
    mode=finalPowerEPS.componentThermalOrbit(out,cfg,a,struct('showSummary',false));
    for label={'nominal','peak'}
        v=mode.(label{1});
        check(max(abs(v.energyResidual_W))<1e-7,[label{1} ' cross-component energy closes.']);
        check(max(abs(v.solar.balanceResidual_W))<1e-7,[label{1} ' two-node thermal energy closes.']);
        check(all(v.solar.electricalPower_W<=v.solar.potentialElectricalPower_W+1e-7),[label{1} ' power extraction stays feasible.']);
        check(~v.solar.structureExchangeEnabled,[label{1} ' unknown structure paths remain disabled.']);
    end
    check(mode.params.downlinkDuration_s==7.55*60 && mode.params.cameraDutyFraction==1,'Uploaded contact and camera assumptions are preserved.');
end
fprintf('PASS: %d analytic, electrical-calibration and integration checks.\n',count);
    function check(ok,message)
        assert(all(ok(:)),message); count=count+1;
    end
    function near(value,expected,tolerance,message)
        check(isequal(size(value),size(expected)) && all(abs(value(:)-expected(:))<=tolerance),message);
    end
end
