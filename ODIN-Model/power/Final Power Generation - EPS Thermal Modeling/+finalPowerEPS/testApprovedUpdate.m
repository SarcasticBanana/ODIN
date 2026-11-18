function results=testApprovedUpdate()
% TESTAPPROVEDUPDATE Independent checks for approved timing/environment additions.
cfg=finalPowerEPS.config(); cfg.h=600; out=finalPowerEPS.propagate(cfg,1);
assert(abs(out.t(end)-out.T)<1e-9,'Exact orbit endpoint missing.');
assert(abs(out.a-(cfg.Re+cfg.h))<1e-9,'Orbit metadata must retain semi-major axis.');
assert(all(diff(out.t)>0),'Event insertion must preserve unique increasing times.');
contacts=[out.eclipseTimes.umbraIn,out.eclipseTimes.umbraOut, ...
    out.eclipseTimes.penumbraIn,out.eclipseTimes.penumbraOut];
for contact=contacts
    assert(min(abs(out.t-contact))<1e-9,'An eclipse contact is absent from the thermal grid.');
end
assert(all(isfinite(out.illum)) && all(out.illum>=0 & out.illum<=1));
longOut=finalPowerEPS.propagate(cfg,10);
assert(all(longOut.illum>=0 & longOut.illum<=1),'Refined disk contacts must not produce negative or >100% illumination.');

p=finalPowerEPS.arrayPower(); radius=cfg.Re+cfg.h;
probe.t=[0 1]; probe.r=[radius -radius;0 0;0 0];
probe.rsun=repmat([cfg.AU;0;0],1,2);
probe.panelNormal_ECI=[-1 1;0 0;0 0];
e=finalPowerEPS.panelEnvironment(probe,cfg,p,0);
exact=(cfg.Re/radius)^2;
assert(max(abs(e.frontEarthViewFactor-exact))<1e-12,'Nadir disk view must match (Re/r)^2.');
assert(all(e.backEarthViewFactor==0),'Opposite face cannot see a nadir disk.');
assert(max(abs(e.frontEarthIR_W_m2-p.frontEmissivity*p.environment.earthInfrared_W_m2*exact))<1e-10);
assert(e.frontAlbedo_W_m2(1)>0 && e.frontAlbedo_W_m2(2)==0,'Only illuminated visible terrain reflects sunlight.');
assert(all(e.backAbsorbedEarthFlux_W_m2==0));
for method={'linear','previous','next','nearest'}
    for query=[0 .17 .5 1]
        q=finalPowerEPS.panelEnvironmentAt(p,e,probe.t,query,method{1});
        for field={'frontEarthViewFactor','backEarthViewFactor','frontAbsorbedEarthFlux_W_m2','backAbsorbedEarthFlux_W_m2'}
            expected=interp1(probe.t,e.(field{1}).',query,method{1}).';
            assert(max(abs(q.(field{1})-expected))<1e-10,'Fast scalar interpolation changed a boundary value.');
        end
    end
end

% A horizon-oriented face exercises partial-disk clipping and quadrature convergence.
probe.panelNormal_ECI=repmat([0;1;0],1,2);
coarse=finalPowerEPS.panelEnvironment(probe,cfg,p,0);
fineP=p; fineP.environment.radialOrder=48; fineP.environment.azimuthCount=192;
fine=finalPowerEPS.panelEnvironment(probe,cfg,fineP,0);
quadratureError=max(abs(coarse.frontAbsorbedEarthFlux_W_m2-fine.frontAbsorbedEarthFlux_W_m2));
assert(quadratureError<0.2,'Earth disk quadrature has not converged sufficiently.');

quiet=struct('showSummary',false,'warnAssumptions',false, ...
    'frontEmissivity',0,'backEmissivity',0,'initialTemperature_C',20,'initialBackTemperature_C',20);
s=finalPowerEPS.SolarArrayThermalModel(struct('irradiance_W_m2',0,'cellTemperature_C',20, ...
    'backTemperature_C',20,'frontAbsorbedEarthFlux_W_m2',0, ...
    'backAbsorbedEarthFlux_W_m2',10),quiet);
assert(abs(s.temperatureRate_K_s)<1e-12 && s.backTemperatureRate_K_s>0, ...
    'Rear Earth heat must enter the rear thermal node.');
assert(abs(s.QheatInput_W-10*s.params.area_m2)<1e-10);
q=quiet; q.frontEmissivity=.85; q.backEmissivity=.85;
s=finalPowerEPS.SolarArrayThermalModel(struct('timeEdges_s',0:10:600,'irradiance_W_m2',0, ...
    'frontEarthViewFactor',.2,'backEarthViewFactor',.4, ...
    'frontAbsorbedEarthFlux_W_m2',20,'backAbsorbedEarthFlux_W_m2',60),q);
assert(max(abs(s.balanceResidual_W))<1e-7,'Earth-inclusive transient must conserve energy.');
assert(max(abs(s.Qout_W-s.QradiationToSpace_W-s.QradiationToEarth_W-s.QtoStructure_W))<1e-9);
assert(max(abs(s.Qradiation_W-s.QfrontRadiation_W-s.QbackRadiation_W))<1e-9);

% Radiation-channel partition must not double physical emitting area.
q.spaceTemperature_K=0; % equal zero-temperature receivers isolate emitting area
base=finalPowerEPS.SolarArrayThermalModel(struct('cellTemperature_C',60,'backTemperature_C',40,'irradiance_W_m2',0),q);
split=finalPowerEPS.SolarArrayThermalModel(struct('cellTemperature_C',60,'backTemperature_C',40, ...
    'irradiance_W_m2',0,'frontEarthViewFactor',.2,'backEarthViewFactor',.4),q);
assert(abs(base.Qout_W-split.Qout_W)<1e-8,'Changing receiver partitions cannot double radiating area.');

% Cached fits must invalidate when calibration changes and restore correctly.
ivInput=struct('time_s',0,'irradiance_W_m2',1353,'cellTemperature_C',28);
ivOptions=struct('mode','xtj_reference','makePlots',false,'showSummary',false,'seriesCells',8);
iv1=finalPowerEPS.solarArrayIV(ivInput,ivOptions);
ivOptions.cellVmp_V=2.40; iv2=finalPowerEPS.solarArrayIV(ivInput,ivOptions);
assert(abs(iv2.Vmp_V-19.2)<1e-7);
ivOptions=rmfield(ivOptions,'cellVmp_V'); iv3=finalPowerEPS.solarArrayIV(ivInput,ivOptions);
assert(abs(iv3.Vmp_V-iv1.Vmp_V)<1e-10 && abs(iv3.Vmp_V-19.248)<1e-7);

results=struct('endpointError_s',out.t(end)-out.T,'horizonFluxRefinement_W_m2',quadratureError, ...
    'thermalResidual_W',max(abs(s.balanceResidual_W)));
fprintf('PASS: timing, analytic Earth geometry, albedo day/night, node assignment, radiation partition, energy conservation and fit-cache invalidation.\n');
disp(results);
end
