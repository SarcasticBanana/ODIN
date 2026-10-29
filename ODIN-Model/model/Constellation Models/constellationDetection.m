%% Constellation Debris Detection Trade Study
%
% Version: 2.5 - MASTER altitude-diameter density and sensor response
% Based on Sean Elliott's V2.3. Last modified: 2026-10-28.
%
% Purpose:
%   Estimate qualified OBJECT PASSAGES per day and over the mission, using
%   a fixed debris environment. A passage is not a unique object or a sensor
%   pixel event. Constellation totals include observations of the same object
%   by different spacecraft. No individual debris orbits are propagated.
%
% Changes relative to V2.3 (updated in V2.5):
%   - Import one broad-range MASTER altitude/diameter density table.
%   - Sum its diameter cells into each configured size bin at each altitude.
%   - Ignore declination dependence; retain V2.3 inclination/RAAN sampling.
%   - Keep the original circular-orbit inclination and relative-RAAN weights.
%   - Interpolate minimum residence time from the selected sensor CSV.
%   - Calculate straight-line FOV chord lengths, rather than qualify a mean
%     crossing time. Numerical altitude layers do not create extra passages.
%   - Include comms downtime, Earth limb clearance, and Sun exclusion.
%   - Report science-size and outside-science-size qualified passages.
%   - Remove unique/repeat estimates and the terminator cosine multiplier.
%
% Intended fidelity / important assumptions:
%   This remains an architecture-level model, NOT a debris propagator or an
%   optical detection simulator. MASTER density is averaged over declination
%   and right ascension. This assumes spatial uniformity on each altitude
%   shell. The independent
%   V2.3 inclination distribution and spatial density are not a jointly
%   fitted phase-space model.
%   Velocities are local T-N approximations at mutual orbit-plane crossings.
%   FOV crossings are straight lines with constant local velocity. Density,
%   sensor response and optional phase angle are evaluated at the midpoint
%   of each observable chord; this approximation should be checked when the
%   response varies rapidly over a crossing. Real pointing is held locally
%   fixed during each crossing. There is no extra angular-rate/velocity cut.
%   Illumination is already represented by the supplied sensor response.
%   No albedo rescaling, terminator factor, or additional shadow penalty is
%   applied. A fixed-phase CSV is not valid at every mission phase angle.
%   MASTER is held fixed over the mission; mission counts are an extrapolation
%   of representative samples, not a forecast of population evolution.
%
% Units:
%   km, seconds, degrees, and objects/km^3 unless a suffix says otherwise.
%   Diameter inputs use cm; MASTER diameter headers use m.
%
% Contributors:
%   Sean Elliott - original V2.3 architecture and inclination sampling.
%
clear;
clc;
close all;
modelFolder = fileparts(mfilename('fullpath'));

%% 1. USER INPUTS

% ---------------- Mission ----------------
Mission.StartDate = datetime(2027,1,1);
Mission.Duration_days = 5*365.25;

% ---------------- Debris Shell ----------------
Shell.MinAltitude_km = 600;
Shell.MaxAltitude_km = 1000;

% ---------------- Spacecraft Orbit ----------------
Orbit.Altitude_km = 600;
Orbit.Type = "SSO";                    % "SSO" or "Circular"
Orbit.Inclination_deg = 97.36;         % Used only for "Circular"
Orbit.NominalLTAN_hr = 18;
Orbit.NumPlanes = 3;
Orbit.SpacecraftPerPlane = 6;
Orbit.BaselineLTANSpread_min = 15;

% ---------------- Sensor Geometry / Response ----------------
Sensor.CSVFile = fullfile(modelFolder,'Sensor_Residence_Time_Trade.csv');
Sensor.FOVShape = "Square";         % "Circular" or "Square"
Sensor.UseCSVFOV = true;
Sensor.FOV_deg = 1.4;                 % FULL angle; used if UseCSVFOV=false
% CSV HFOV_deg is treated as FULL width. For Square this is each side's
% angular width, not the diagonal. For Circular it is full cone angle.
Sensor.ExpectedAlbedo = 0.20;         % Must agree with CSV, not a multiplier
Sensor.UsePhaseAngleCSV = false;      % Enable only with a phase-angle sweep
Sensor.FixedPhaseAngle_deg = 30;      % Response slice when toggle is false
Sensor.ReportPhaseDiagnostic = false;
% Range is taken from CSV coverage, including its points below 50 km.
% Range=0 is a mathematical endpoint only. No range/diameter/phase extrapolation.

% ---------------- Pointing ----------------
% "Zenith" -> +R, "AlongTrack" -> +T, "Normal+" -> +N,
% "Normal-" -> -N, or "CustomRTN" -> normalized [R T N].
Pointing.Mode = "CustomRTN";
Pointing.CustomRTN = [sind(66) 0 cosd(66)];
% 66 deg FROM orbit normal TOWARD zenith. RTN direction is not phase angle.
% Choose the N sign opposite the Sun at each representative geometry. The R
% and T components are unchanged. This does not guarantee Sun clearance.
Pointing.ChooseAntiSunNormal = false;

% ---------------- Debris Population ----------------
% Use the broad altitude-diameter export, not separate size-filtered cases.
% Filename is only a selector; column meanings and diameter bounds come from
% the file header. The 1-7 cm export was checked against the earlier 1-6 cm
% export: retained science bins agree, and the former 4.95 cm spike disappears.
% This is a practical input choice, not a validation of MASTER's internals.
Debris.DensityModel = "MASTER";        % "MASTER" or "Demonstration"
Debris.MASTERFolder = fullfile(modelFolder,'MASTER_Data');
Debris.MASTERFileName = "ODIN_2024_D010_070_t.__1";
% The folder and its subfolders are searched automatically for this filename.
% Set MASTERFileName="" to discover an altitude-diameter table by its header;
% auto-discovery requires exactly one suitable table, never a guessed latest
% file. Old altitude-declination tables are ignored in auto-discovery mode.
Debris.SizeBins_cm = [1 1.5; 1.5 2.5; 2.5 3.1; 3.1 5];
Debris.RepresentativeDiameter_cm = [1.25 2.5 3 5];
Debris.ScienceSizeRange_cm = [1 3.1];
Debris.MaxUsableDiameter_cm = 5;
% Data above this explicit limit are not counted. For a larger modeled size
% range, review the export's upper-bound behavior and sensor CSV coverage,
% then deliberately change the limit and SizeBins_cm together.
% Density belongs to the ENTIRE bin. Only the sensor response uses its
% representative diameter. Diameter edges must coincide with MASTER cell
% edges: partial cells are rejected rather than split by an assumed size law.
Debris.DemoDensity_per_km3 = [2e-7 1e-7 4e-8 3e-8]; % Synthetic, NOT predictions
Debris.ExpectedReferenceEpoch = "";   % Optional exact header, e.g. "2024/08"
% The selected population snapshot is held fixed throughout the mission.
% There is no declination weighting and no size-bin renormalization.

% Same inclination and relative-RAAN sampling as V2.3.
% Uniform in inclination angle is NOT isotropic in orbit-plane orientation.
Debris.InclinationModel = "Uniform";  % "Uniform" or "Function"
Debris.InclinationFunction = @(i_deg) ones(size(i_deg));
Debris.NumInclinationSamples = 100;
Debris.NumRAANSamples = 8;

% ---------------- Operations / Exclusion Angles ----------------
Comms.Downlink_min_per_day = 120;     % Total downlink time across all sessions
Comms.Transition_min_per_session = 10;% Total overhead per session, not per end
Comms.Sessions_per_day = 1;


% Availability is applied ONCE after geometric detection rates. Until a
% schedule is supplied, downtime is assumed independent of orbital geometry.
Exclusion.SunClearance_deg = 5;       % Extra clearance from solar disk to FOV
Exclusion.EarthLimbClearance_deg = 0.05;
% Earth masks rays against the apparent Earth disk plus margin. This also
% excludes foreground objects against the Earth, not just occulted objects.
% Sun blocks the ENTIRE FOV when the exclusion disk can touch it. For Square,
% a circumscribed cone is used (conservative at the corners).

% ---------------- Architecture Trades ----------------
Trade.RunArchitectureTrade = false;  % Enable once baseline inputs are checked
Trade.NumPlanes = 1:6;
Trade.SpacecraftPerPlane = 1:9;
Trade.ArchitectureLTANSpread_min = 60;
Trade.RunLTANTrade = false;
Trade.LTANSpread_min = 0:10:180;
Trade.RunAltitudeTrade = false;
Trade.SpacecraftAltitude_km = 570:10:700;
% Optional trades retain V2.3's intent, but no unique/repeat-object plots.

% ---------------- Numerical Settings ----------------
Numerics.VelocityAltitudeStep_km = 50;
% Velocity is computed at each altitude interval's center using the
% circular-orbit formula. Whole FOV chords are then assigned by their midpoint
% altitude. Internal altitude boundaries do NOT truncate the residence time.
Numerics.ProjectedAlongSamples = 40;
Numerics.ProjectedAcrossSamples = 12;
Numerics.OrbitalPositionSamples = 12;
Numerics.MissionDateSamples = 24;
Numerics.RunGeometrySelfChecks = true; % Small analytic checks before the main run
% Midpoint samples span the full mission duration and orbital argument of
% latitude. Increase date/orbit/area/orientation resolution for convergence.
% These sample counts are a starting setting, not a demonstrated convergence.

% ---------------- Plot / Output Controls ----------------
Plot.MakeFigures = true;
Plot.SaveFigures = false;
Plot.OutputFolder = fullfile(modelFolder,'trade_study_figures');
Plot.ExportResolution_dpi = 300;
Output.SaveTables = false;
Output.Folder = fullfile(modelFolder,'model_results');
Output.ReportPoissonProbability = false;
% Optional P(at least one) assumes a stationary Poisson passage process;
% this is NOT probability of detecting a distinct object in the population.

%% 2. CONSTANTS AND INPUT CHECKS

Earth.Radius_km = 6371;
Earth.Mu_km3_s2 = 398600;
Earth.J2 = 1.08262668e-3;
secondsPerDay = 86400;
tropicalYear_days = 365.2422;
assert(Mission.Duration_days>0 && isfinite(Mission.Duration_days));
assert(Shell.MaxAltitude_km>Shell.MinAltitude_km && Shell.MinAltitude_km>0);
assert(Orbit.NumPlanes>=1 && mod(Orbit.NumPlanes,1)==0);
assert(Orbit.SpacecraftPerPlane>=1 && mod(Orbit.SpacecraftPerPlane,1)==0);
assert(all([Comms.Downlink_min_per_day Comms.Transition_min_per_session ...
    Comms.Sessions_per_day]>=0),'Comms values must be nonnegative.');
Comms.TotalDowntime_min_per_day = Comms.Downlink_min_per_day + ...
    Comms.Transition_min_per_session*Comms.Sessions_per_day;
assert(Comms.TotalDowntime_min_per_day<=1440,'Downtime exceeds one day.');
Comms.Availability = 1-Comms.TotalDowntime_min_per_day/1440;
assert(all([Exclusion.SunClearance_deg Exclusion.EarthLimbClearance_deg]>=0));
validateattributes(Debris.SizeBins_cm,{'numeric'},{'2d','ncols',2,'finite','positive'});
nBins = size(Debris.SizeBins_cm,1);
assert(numel(Debris.RepresentativeDiameter_cm)==nBins);
assert(all(Debris.SizeBins_cm(:,2)>Debris.SizeBins_cm(:,1)));
assert(all(Debris.SizeBins_cm(2:end,1)>=Debris.SizeBins_cm(1:end-1,2)), ...
    'Size bins must be sorted and disjoint.');
for k = 1:nBins
    assert(Debris.RepresentativeDiameter_cm(k)>=Debris.SizeBins_cm(k,1) && ...
        Debris.RepresentativeDiameter_cm(k)<=Debris.SizeBins_cm(k,2));
end
science = Debris.SizeBins_cm(:,1)>=Debris.ScienceSizeRange_cm(1) & ...
    Debris.SizeBins_cm(:,2)<=Debris.ScienceSizeRange_cm(2);
straddle = Debris.SizeBins_cm(:,1)<Debris.ScienceSizeRange_cm(2) & ...
    Debris.SizeBins_cm(:,2)>Debris.ScienceSizeRange_cm(1) & ~science;
assert(~any(straddle),'Split bins at the science-size boundaries.');
assert(any(science),'No bin is inside the science-size interval.');
assert(all(Debris.SizeBins_cm(:,2)<=Debris.MaxUsableDiameter_cm), ...
    'A configured bin exceeds MaxUsableDiameter_cm.');
validateattributes(Numerics.VelocityAltitudeStep_km,{'numeric'},{'scalar','positive','finite'});
counts = [Numerics.ProjectedAlongSamples Numerics.ProjectedAcrossSamples ...
    Numerics.OrbitalPositionSamples Numerics.MissionDateSamples ...
    Debris.NumInclinationSamples Debris.NumRAANSamples];
assert(all(counts>=2 & mod(counts,1)==0),'Sample counts must be integers >=2.');

%% 3. LOAD SENSOR RESPONSE AND MASTER ENVIRONMENT

if Numerics.RunGeometrySelfChecks, runGeometrySelfChecks; end
[Response,Sensor] = readSensorResponse(Sensor,Debris);
Environment = readMASTEREnvironment(Debris,Shell,Earth);
PopulationSummary = Environment.Summary;
disp(PopulationSummary);
Orbit = buildOrbitAtAltitude(Orbit.Altitude_km,Orbit,Earth, ...
    tropicalYear_days,secondsPerDay);
if strcmpi(Debris.DensityModel,'Demonstration')
    runLabel = "SYNTHETIC DEMONSTRATION - NOT A MISSION PREDICTION";
    warning('%s',runLabel);
else
    runLabel = "MASTER simulation " + string(Environment.SimulationInterval) + ...
        " (reference " + Environment.ReferenceEpoch + ")";
    fprintf('MASTER altitude-diameter table loaded; grid and coverage checks passed.\n');
    fprintf('No independent population or separate-run consistency validation is claimed.\n');
end
if ~Sensor.UsePhaseAngleCSV
    warning(['Fixed %.1f deg phase response is used at every orbital sample. ' ...
        'It is not validated for all illumination conditions.'],Sensor.FixedPhaseAngle_deg);
end
fprintf('\n%s\n',runLabel);
fprintf('FOV: %s, full width %.6f deg; CSV range %.3f to %.3f km.\n', ...
    Sensor.FOVShape,Sensor.FOV_deg,Response.Range_km(1),Response.Range_km(end));
fprintf('Downtime: %.1f min/day; availability %.6f.\n', ...
    Comms.TotalDowntime_min_per_day,Comms.Availability);

%% 4. BASELINE ENCOUNTER GEOMETRY AND DETECTION RATE

% A projected-area quadrature samples straight paths through the FOV.
% Each accepted path contributes rho * v_rel * dA * orientation_weight.
% An observable chord is clipped only by the physical shell, true CSV range,
% FOV, and Earth exclusion. Adjacent numerical density layers are not separate
% detections. Disconnected observable intervals count as separate passages.
% Density depends on chord-midpoint altitude and diameter bin only.
% Representative spacecraft positions still determine Sun exclusion and
% optional phase angle; they do not introduce declination weighting.
Geometry = buildEncounterGeometry(Orbit,Shell,Sensor,Pointing,Debris, ...
    Exclusion,Numerics,Earth);
offsets_min = distributePlanesInLTAN(Orbit.NumPlanes,Orbit.BaselineLTANSpread_min);
[planeRate_s,Diagnostics] = evaluatePlanes(offsets_min,Geometry,Orbit, ...
    Mission,Sensor,Response,Environment,Debris,Pointing,Exclusion,Numerics,Earth);
planeRate_s = planeRate_s*Comms.Availability;
spacecraftDay = mean(planeRate_s,1)'*secondsPerDay;
constellationDay = sum(planeRate_s,1)'*Orbit.SpacecraftPerPlane*secondsPerDay;

%% 5. BASELINE RESULTS

% Per-spacecraft results are averaged across planes. PlaneResults preserves
% each plane's individual daily rate when LTAN geometry differs.
Results = table(Debris.SizeBins_cm(:,1),Debris.SizeBins_cm(:,2), ...
    Debris.RepresentativeDiameter_cm(:),science,spacecraftDay,constellationDay, ...
    spacecraftDay*Mission.Duration_days,constellationDay*Mission.Duration_days, ...
    'VariableNames',{'LowerDiameter_cm','UpperDiameter_cm','Representative_cm', ...
    'ScienceBin','PassagesPerSpacecraftPerDay','PassagesConstellationPerDay', ...
    'PassagesPerSpacecraftMission','PassagesConstellationMission'});
category = ["Science sizes";"Outside science sizes";"All modeled sizes"];
groupMasks = [science ~science true(nBins,1)];
scGroup = groupMasks'*spacecraftDay;
conGroup = groupMasks'*constellationDay;
Summary = table(category,scGroup,conGroup,scGroup*Mission.Duration_days, ...
    conGroup*Mission.Duration_days,'VariableNames',{'Category', ...
    'PassagesPerSpacecraftPerDay','PassagesConstellationPerDay', ...
    'PassagesPerSpacecraftMission','PassagesConstellationMission'});
if Output.ReportPoissonProbability
    Summary.P_AtLeastOne_ConstellationDay = -expm1(-conGroup);
end
PlaneResults = array2table([offsets_min(:),planeRate_s*secondsPerDay], ...
    'VariableNames',[{'LTANOffset_min'},cellstr(compose('Bin%d_PerSpacecraftPerDay',1:nBins))]);
Results.EnvironmentMode = repmat(string(Debris.DensityModel),nBins,1);
Summary.EnvironmentMode = repmat(string(Debris.DensityModel),height(Summary),1);
PlaneResults.EnvironmentMode = repmat(string(Debris.DensityModel),height(PlaneResults),1);
disp(Results);
disp(Summary);
fprintf('Counts are qualified passages, not unique debris or sensor pixel events.\n');
fprintf('Total includes ONLY %.3g-%.3g cm covered by the configured bins.\n', ...
    Debris.SizeBins_cm(1,1),Debris.SizeBins_cm(end,2));
fprintf('Mission totals use a fixed environment and sampled mean operating geometry.\n');
if Sensor.ReportPhaseDiagnostic
    fprintf('Sampled chord-midpoint phase span: %.2f to %.2f deg.\n', ...
        Diagnostics.MinPhase_deg,Diagnostics.MaxPhase_deg);
end

%% 6. OPTIONAL ARCHITECTURE, LTAN, AND ALTITUDE TRADES

% Spacecraft-per-plane scaling is linear because these are passages, not
% deduplicated objects. Plane-count/LTAN changes recompute sampled geometry.
ArchitectureResults = table;
LTANResults = table;
AltitudeResults = table;
if Trade.RunArchitectureTrade
    rows = [];
    for p = Trade.NumPlanes
        offsets = distributePlanesInLTAN(p,Trade.ArchitectureLTANSpread_min);
        rates = evaluatePlanes(offsets,Geometry,Orbit,Mission,Sensor,Response, ...
            Environment,Debris,Pointing,Exclusion,Numerics,Earth);
        for n = Trade.SpacecraftPerPlane
            daily = sum(rates,1)*n*secondsPerDay*Comms.Availability;
            rows(end+1,:) = [p n sum(daily(science)) sum(daily(~science)) sum(daily)]; %#ok<SAGROW>
        end
    end
    ArchitectureResults = array2table(rows,'VariableNames', ...
        {'Planes','SpacecraftPerPlane','SciencePerDay','OutsidePerDay','TotalPerDay'});
end
if Trade.RunLTANTrade
    rows = [];
    for spread = Trade.LTANSpread_min
        offsets = distributePlanesInLTAN(Orbit.NumPlanes,spread);
        rates = evaluatePlanes(offsets,Geometry,Orbit,Mission,Sensor,Response, ...
            Environment,Debris,Pointing,Exclusion,Numerics,Earth);
        daily = sum(rates,1)*Orbit.SpacecraftPerPlane*secondsPerDay*Comms.Availability;
        rows(end+1,:) = [spread sum(daily(science)) sum(daily(~science)) sum(daily)]; %#ok<SAGROW>
    end
    LTANResults = array2table(rows,'VariableNames', ...
        {'LTANSpread_min','SciencePerDay','OutsidePerDay','TotalPerDay'});
end
if Trade.RunAltitudeTrade
    rows = [];
    for h = Trade.SpacecraftAltitude_km
        caseOrbit = buildOrbitAtAltitude(h,Orbit,Earth,tropicalYear_days,secondsPerDay);
        caseGeometry = buildEncounterGeometry(caseOrbit,Shell,Sensor,Pointing, ...
            Debris,Exclusion,Numerics,Earth);
        rates = evaluatePlanes(offsets_min,caseGeometry,caseOrbit,Mission,Sensor, ...
            Response,Environment,Debris,Pointing,Exclusion,Numerics,Earth);
        daily = sum(rates,1)*Orbit.SpacecraftPerPlane*secondsPerDay*Comms.Availability;
        rows(end+1,:) = [h sum(daily(science)) sum(daily(~science)) sum(daily)]; %#ok<SAGROW>
    end
    AltitudeResults = array2table(rows,'VariableNames', ...
        {'Altitude_km','SciencePerDay','OutsidePerDay','TotalPerDay'});
end

%% 7. FIGURES AND SAVED OUTPUTS

if Plot.MakeFigures
    f = figure('Color','w','Name','Qualified object passages');
    bar([spacecraftDay constellationDay]);
    xticks(1:nBins);
    xticklabels(compose('%.3g-%.3g cm',Debris.SizeBins_cm(:,1),Debris.SizeBins_cm(:,2)));
    ylabel('Qualified object passages / day');
    legend('Mean spacecraft','Whole constellation','Location','best');
    title({'Qualified passages by size',char(runLabel)});
    grid on;
    if Plot.SaveFigures
        if ~isfolder(Plot.OutputFolder), mkdir(Plot.OutputFolder); end
        exportgraphics(f,fullfile(Plot.OutputFolder,'qualified_passages.png'), ...
            'Resolution',Plot.ExportResolution_dpi);
    end
    plotTradeResults(ArchitectureResults,LTANResults,AltitudeResults,runLabel);
end
if Output.SaveTables
    if ~isfolder(Output.Folder), mkdir(Output.Folder); end
    % Separate run folder prevents a demonstration from overwriting real results.
    runFolder = fullfile(Output.Folder,char(datetime('now','Format','yyyyMMdd_HHmmss_SSS')));
    mkdir(runFolder);
    writetable(Results,fullfile(runFolder,'size_bin_passages.csv'));
    writetable(PopulationSummary,fullfile(runFolder,'imported_population_summary.csv'));
    writetable(Summary,fullfile(runFolder,'summary_passages.csv'));
    writetable(PlaneResults,fullfile(runFolder,'plane_passages.csv'));
    if ~isempty(ArchitectureResults), writetable(ArchitectureResults,fullfile(runFolder,'architecture_trade.csv')); end
    if ~isempty(LTANResults), writetable(LTANResults,fullfile(runFolder,'LTAN_trade.csv')); end
    if ~isempty(AltitudeResults), writetable(AltitudeResults,fullfile(runFolder,'altitude_trade.csv')); end
    save(fullfile(runFolder,'run_inputs_and_results.mat'),'Mission','Shell','Orbit', ...
        'Sensor','Pointing','Debris','Comms','Exclusion','Numerics','Trade', ...
        'Results','Summary','PlaneResults','PopulationSummary','Environment','Diagnostics','runLabel');
    fid = fopen(fullfile(runFolder,'RUN_STATUS.txt'),'w');
    assert(fid>=0,'Cannot write run status.');
    fprintf(fid,'%s\nEnvironment mode: %s\nSensor CSV: %s\n',runLabel, ...
        Debris.DensityModel,Sensor.CSVFile);
    fprintf(fid,'Fixed environment; qualified passages, not unique objects.\n');
    fprintf(fid,'Altitude-diameter density; spatially uniform over each shell.\n');
    fprintf(fid,'Only configured bins through %.6g cm are counted.\n',Debris.MaxUsableDiameter_cm);
    fprintf(fid,'Phase interpolation enabled: %d\n',Sensor.UsePhaseAngleCSV);
    if isfield(Environment,'SourceFiles')
        fprintf(fid,'MASTER source: %s\n',Environment.SourceFiles{:});
    end
    fclose(fid);
    fprintf('Saved results: %s\n',runFolder);
end

%% 8. LOCAL FUNCTIONS

function [Response,Sensor] = readSensorResponse(Sensor,Debris)
% Read a complete diameter/range grid (and optionally phase grid).
% MaxRelativeVel_kms and MaxAngularRate_degs are intentionally NOT thresholds.
% They may be infinite at zero range without invalidating the residence table.
    assert(isfile(Sensor.CSVFile),'Sensor CSV not found: %s',Sensor.CSVFile);
    T = readtable(Sensor.CSVFile);
    required = {'DebrisSize_cm','SlantRange_km','MinResidenceTime_s','PhaseAngle_deg'};
    assert(all(ismember(required,T.Properties.VariableNames)), ...
        'CSV requires DebrisSize_cm, SlantRange_km, MinResidenceTime_s, PhaseAngle_deg.');
    assert(all(isfinite(T{:,required}),'all'),'Nonfinite sensor response entries.');
    assert(all(T.DebrisSize_cm>0 & T.SlantRange_km>=0 & T.MinResidenceTime_s>=0));
    assert(all(T.PhaseAngle_deg>=0 & T.PhaseAngle_deg<=180));
    if ismember('Albedo',T.Properties.VariableNames)
        assert(all(abs(T.Albedo-Sensor.ExpectedAlbedo)<1e-8), ...
            'CSV albedo differs from ExpectedAlbedo. Obtain a matching response; do not rescale it.');
    else
        warning('TBD: CSV does not document its albedo.');
    end
    if Sensor.UseCSVFOV
        assert(ismember('HFOV_deg',T.Properties.VariableNames),'CSV needs HFOV_deg or disable UseCSVFOV.');
        assert(all(isfinite(T.HFOV_deg)) && max(T.HFOV_deg)-min(T.HFOV_deg)<1e-8);
        Sensor.FOV_deg = T.HFOV_deg(1);
    end
    assert(Sensor.FOV_deg>0 && Sensor.FOV_deg<90);
    assert(any(strcmpi(Sensor.FOVShape,{'Circular','Square'})),'Unknown FOV shape.');
    if ~Sensor.UseCSVFOV && ismember('HFOV_deg',T.Properties.VariableNames)
        warning('Manual FOV override: confirm the sensor response remains valid for this FOV.');
    end
    Response.Range_km = unique(T.SlantRange_km);
    Response.Diameter_cm = unique(T.DebrisSize_cm);
    Response.Phase_deg = unique(T.PhaseAngle_deg);
    assert(numel(Response.Range_km)>=2 && numel(Response.Diameter_cm)>=2);
    assert(all(Debris.RepresentativeDiameter_cm>=min(Response.Diameter_cm) & ...
        Debris.RepresentativeDiameter_cm<=max(Response.Diameter_cm)), ...
        'Representative diameters exceed sensor CSV coverage.');
    [~,ir] = ismember(T.SlantRange_km,Response.Range_km);
    [~,id] = ismember(T.DebrisSize_cm,Response.Diameter_cm);
    [~,ip] = ismember(T.PhaseAngle_deg,Response.Phase_deg);
    sz = [numel(Response.Range_km) numel(Response.Diameter_cm) numel(Response.Phase_deg)];
    idx = sub2ind(sz,ir,id,ip);
    assert(numel(unique(idx))==height(T),'Duplicate diameter/range/phase CSV rows.');
    assert(height(T)==prod(sz),'CSV must have a complete rectangular range/diameter/phase grid.');
    values = nan(sz); values(idx) = T.MinResidenceTime_s;
    if Sensor.UsePhaseAngleCSV
        assert(numel(Response.Phase_deg)>=2,'Phase toggle needs a CSV with multiple phase angles.');
        Response.Time = griddedInterpolant({Response.Range_km,Response.Diameter_cm, ...
            Response.Phase_deg},values,'linear','none');
    else
        assert(Sensor.FixedPhaseAngle_deg>=min(Response.Phase_deg) && ...
            Sensor.FixedPhaseAngle_deg<=max(Response.Phase_deg),'Fixed phase is outside CSV coverage.');
        if numel(Response.Phase_deg)==1
            slice = values(:,:,1);
        else
            temp = griddedInterpolant({Response.Range_km,Response.Diameter_cm, ...
                Response.Phase_deg},values,'linear','none');
            [r,d] = ndgrid(Response.Range_km,Response.Diameter_cm);
            slice = temp(r,d,Sensor.FixedPhaseAngle_deg*ones(size(r)));
        end
        Response.Time = griddedInterpolant({Response.Range_km,Response.Diameter_cm},slice,'linear','none');
    end
    Sensor.Range_km = max(Response.Range_km);
    Sensor.MinRange_km = min(Response.Range_km);
end

function Environment = readMASTEREnvironment(Debris,Shell,Earth)
% Load ONE broad-range differential density table with altitude and diameter
% axes. MASTER reports density per diameter cell in objects/km^3, not a density
% per metre of diameter. Sum cells; do NOT multiply by their width again.
% Spatial density is averaged over the exported angular region, required
% here to cover all declinations and right ascensions. No declination model
% is inferred from this table.
    Environment.Mode = Debris.DensityModel;
    nBins = size(Debris.SizeBins_cm,1);
    if strcmpi(Debris.DensityModel,'Demonstration')
        assert(numel(Debris.DemoDensity_per_km3)==nBins);
        assert(all(isfinite(Debris.DemoDensity_per_km3) & Debris.DemoDensity_per_km3>=0));
        Environment.DemoDensity = Debris.DemoDensity_per_km3;
        Environment.ReferenceEpoch = "SYNTHETIC";
        values = Debris.DemoDensity_per_km3(:);
        Environment.Summary = table(Debris.SizeBins_cm(:,1),Debris.SizeBins_cm(:,2), ...
            values,values,values,values,zeros(nBins,1), ...
            'VariableNames',{'LowerDiameter_cm','UpperDiameter_cm', ...
            'MinDensity_per_km3','MedianDensity_per_km3','MaxDensity_per_km3', ...
            'ShellMeanDensity_per_km3','ZeroAltitudeCells'});
        return;
    end
    assert(strcmpi(Debris.DensityModel,'MASTER'),'Unknown density model.');
    assert(isfolder(Debris.MASTERFolder),'MASTER input folder does not exist: %s',Debris.MASTERFolder);
    if strlength(Debris.MASTERFileName)>0
        files = dir(fullfile(Debris.MASTERFolder,'**',char(Debris.MASTERFileName)));
        files = files(~[files.isdir]);
        assert(numel(files)==1, ...
            'Expected one copy of %s beneath %s; found %d.', ...
            Debris.MASTERFileName,Debris.MASTERFolder,numel(files));
        path = fullfile(files(1).folder,files(1).name);
        c = parseMASTERAltitudeDiameter(fileread(path),path);
    else
        files = dir(fullfile(Debris.MASTERFolder,'**','*'));
        candidates = {};
        for k = 1:numel(files)
            if files(k).isdir || isempty(regexp(files(k).name,'\.__\d+$','once')), continue; end
            path = fullfile(files(k).folder,files(k).name);
            raw = fileread(path);
            x = regexp(raw,'x Data:\s*([^\r\n]+)','tokens','once');
            y = regexp(raw,'y Data:\s*([^\r\n]+)','tokens','once');
            if isempty(x) || isempty(y), continue; end
            axes = {strtrim(x{1}),strtrim(y{1})};
            if ~all(ismember({'S.D. Altitude [km]','S.D. Diameter [m]'},axes)), continue; end
            candidate = parseMASTERAltitudeDiameter(raw,path);
            if candidate.DiameterEdges_cm(1)<=min(Debris.SizeBins_cm(:,1))+1e-7 && ...
                    candidate.DiameterEdges_cm(end)>=max(Debris.SizeBins_cm(:,2))-1e-7 && ...
                    candidate.AltitudeEdges_km(1)<=Shell.MinAltitude_km && ...
                    candidate.AltitudeEdges_km(end)>=Shell.MaxAltitude_km
                candidates{end+1} = candidate; %#ok<AGROW>
            end
        end
        assert(numel(candidates)==1, ...
            ['Auto-discovery found %d suitable altitude-diameter tables. ' ...
            'Set MASTERFileName to the intended broad export.'],numel(candidates));
        c = candidates{1};
    end
    if strlength(Debris.ExpectedReferenceEpoch)>0
        assert(strcmp(c.Epoch,Debris.ExpectedReferenceEpoch),'Unexpected MASTER reference epoch.');
    end
    assert(c.AltitudeEdges_km(1)<=Shell.MinAltitude_km && ...
        c.AltitudeEdges_km(end)>=Shell.MaxAltitude_km,'MASTER grid does not cover the modeled shell.');
    assert(min(Debris.SizeBins_cm(:,1))>=c.SelectedDiameter_cm(1)-1e-7 && ...
        max(Debris.SizeBins_cm(:,2))<=c.SelectedDiameter_cm(2)+1e-7, ...
        'Configured diameter bins exceed the selected MASTER population bounds.');
    density = zeros(numel(c.Altitude_km),nBins);
    cellsPerBin = zeros(nBins,1);
    for b = 1:nBins
        lower = find(abs(c.DiameterEdges_cm-Debris.SizeBins_cm(b,1))<1e-7);
        upper = find(abs(c.DiameterEdges_cm-Debris.SizeBins_cm(b,2))<1e-7);
        assert(numel(lower)==1 && numel(upper)==1 && upper>lower, ...
            ['Bin %.6g-%.6g cm must align with MASTER diameter-cell edges. ' ...
            'No extrapolation or partial-cell allocation is performed.'],Debris.SizeBins_cm(b,:));
        density(:,b) = sum(c.CellDensity_per_km3(:,lower:upper-1),2);
        cellsPerBin(b) = upper-lower;
    end
    Environment.ReferenceEpoch = string(c.Epoch);
    Environment.SimulationInterval = c.Interval;
    Environment.SourceFiles = {c.Path};
    Environment.SelectedDiameter_cm = c.SelectedDiameter_cm;
    Environment.Altitude_km = c.Altitude_km;
    Environment.AltitudeEdges_km = c.AltitudeEdges_km;
    Environment.Density_per_km3 = density;
    Environment.DiameterCellsPerBin = cellsPerBin;
    % Report density diagnostics over the modeled shell. Shell mean is only
    % a diagnostic: encounter calculations query the altitude cells directly.
    lo = max(c.AltitudeEdges_km(1:end-1),Shell.MinAltitude_km);
    hi = min(c.AltitudeEdges_km(2:end),Shell.MaxAltitude_km);
    active = hi>lo;
    weights = (Earth.Radius_km+hi(active)).^3-(Earth.Radius_km+lo(active)).^3;
    weights = weights/sum(weights);
    values = density(active,:);
    Environment.Summary = table(Debris.SizeBins_cm(:,1),Debris.SizeBins_cm(:,2), ...
        min(values,[],1)',median(values,1)',max(values,[],1)', ...
        (weights'*values)',sum(values==0,1)',cellsPerBin, ...
        'VariableNames',{'LowerDiameter_cm','UpperDiameter_cm', ...
        'MinDensity_per_km3','MedianDensity_per_km3','MaxDensity_per_km3', ...
        'ShellMeanDensity_per_km3','ZeroAltitudeCells','DiameterCellsSummed'});
    fprintf('\nImported MASTER file: %s\n',c.Path);
    fprintf('Population selection: %.6g-%.6g cm; only configured bins are counted.\n',c.SelectedDiameter_cm);
    if any(all(values==0,1))
        warning('An entire configured diameter bin has zero imported density. See PopulationSummary.');
    end
end

function c = parseMASTERAltitudeDiameter(raw,path)
% Use headers, not filename suffixes, to identify axes. Both axis orders work.
% Cell centers and spacings define edges; declared class count/max labels can
% differ slightly from the actual printed grid, so numeric rows are checked.
    x = regexp(raw,'x Data:\s*([^\r\n]+)','tokens','once');
    y = regexp(raw,'y Data:\s*([^\r\n]+)','tokens','once');
    assert(~isempty(x) && ~isempty(y),'Missing spectrum-axis headers: %s',path);
    axes = {strtrim(x{1}),strtrim(y{1})};
    altitudeColumn = find(strcmp(axes,'S.D. Altitude [km]'));
    diameterColumn = find(strcmp(axes,'S.D. Diameter [m]'));
    assert(numel(altitudeColumn)==1 && numel(diameterColumn)==1, ...
        'Expected altitude-diameter density, not altitude-declination: %s',path);
    assert(~isempty(regexp(raw,'z Data:\s*Spatial Density \[1/km\^3\]','once')), ...
        'Expected spatial-density units [1/km^3]: %s',path);
    assert(~isempty(regexp(raw,'Spectrum:\s*differential','once')), ...
        'Use differential, not cumulative, MASTER output: %s',path);
    number = '([-+]?\d*\.?\d+(?:[EeDd][-+]?\d+)?)';
    lower = regexp(raw,['#\s*' number '[^\r\n]*?: lower diameter threshold'],'tokens','once');
    upper = regexp(raw,['#\s*' number '[^\r\n]*?: upper diameter threshold'],'tokens','once');
    epoch = regexp(raw,'Most recent reference Epoch:\s*([^\r\n]+)','tokens','once');
    interval = regexp(raw,'simulation between\s*([^\r\n]+)','tokens','once');
    assert(~isempty(lower) && ~isempty(upper) && ~isempty(epoch) && ~isempty(interval), ...
        'Missing MASTER population metadata: %s',path);
    c.SelectedDiameter_cm = 100*[str2double(regexprep(lower{1},'[Dd]','E')) ...
        str2double(regexprep(upper{1},'[Dd]','E'))];
    c.Epoch = strtrim(epoch{1}); c.Interval = strtrim(interval{1}); c.Path = path;
    region = regexp(raw,[number '\s*-\s*' number '\s*deg declination and\s*' ...
        number '\s*-\s*' number '\s*deg right ascension'],'tokens','once');
    assert(~isempty(region),'Missing spatial-region bounds in MASTER header: %s',path);
    region = cellfun(@str2double,region);
    assert(max(abs(region-[-90 90 -180 180]))<1e-6, ...
        'Altitude-only model requires full declination and right-ascension coverage.');
    lines = splitlines(string(raw));
    data = [];
    for k = 1:numel(lines)
        line = strtrim(lines(k));
        if strlength(line)==0 || startsWith(line,'#'), continue; end
        values = sscanf(regexprep(char(line),'[Dd]','E'),'%f')';
        assert(numel(values)==20,'Expected two axes and 18 source-density columns: %s',path);
        data(end+1,:) = values; %#ok<AGROW>
    end
    assert(~isempty(data) && all(isfinite(data),'all'),'Missing/nonfinite MASTER numeric data.');
    assert(all(data(:,3:end)>=0,'all'),'Negative MASTER density.');
    assert(all(abs(data(:,20)-data(:,17))<=max(1e-30,abs(data(:,20))*1e-4)), ...
        'Export Human-made Population only; no meteoroid/stream contribution.');
    altitude = data(:,altitudeColumn);
    diameter_cm = 100*data(:,diameterColumn);
    c.Altitude_km = unique(altitude);
    c.Diameter_cm = unique(diameter_cm);
    assert(numel(c.Altitude_km)>=2 && numel(c.Diameter_cm)>=2,'Need at least two cells per axis.');
    dh = diff(c.Altitude_km); dd = diff(c.Diameter_cm);
    assert(max(abs(dh-dh(1)))<1e-6 && max(abs(dd-dd(1)))<1e-7, ...
        'This importer requires uniform linear altitude and diameter grids.');
    c.AltitudeEdges_km = [c.Altitude_km-dh(1)/2;c.Altitude_km(end)+dh(1)/2];
    c.DiameterEdges_cm = [c.Diameter_cm-dd(1)/2;c.Diameter_cm(end)+dd(1)/2];
    [~,ih] = ismember(altitude,c.Altitude_km);
    [~,id] = ismember(diameter_cm,c.Diameter_cm);
    shape = [numel(c.Altitude_km),numel(c.Diameter_cm)];
    indices = sub2ind(shape,ih,id);
    assert(numel(unique(indices))==numel(indices) && numel(indices)==prod(shape), ...
        'Duplicate or missing altitude-diameter cells: %s',path);
    c.CellDensity_per_km3 = zeros(shape);
    c.CellDensity_per_km3(indices) = data(:,17);
end

function density = queryDensity(Environment,bin,altitude)
% Piecewise constant cell-average density. No declination weighting and no
% extrapolation outside the imported altitude coverage.
    if strcmpi(Environment.Mode,'Demonstration')
        density = Environment.DemoDensity(bin)*ones(size(altitude)); return;
    end
    ih = discretize(altitude,Environment.AltitudeEdges_km);
    assert(all(isfinite(ih)),'Requested altitude is outside MASTER grid.');
    density = Environment.Density_per_km3(ih,bin);
end

function OrbitCase = buildOrbitAtAltitude( ...
    altitude_km,OrbitTemplate,Earth,tropicalYear_days,secondsPerDay)
% Recompute circular-orbit properties for one spacecraft altitude.
%
% For an SSO case, inclination is recalculated from the first-order J2
% precession equation so every altitude case remains Sun-synchronous.

    OrbitCase = OrbitTemplate;

    OrbitCase.Altitude_km = altitude_km;
    OrbitCase.Radius_km = Earth.Radius_km + altitude_km;

    OrbitCase.MeanMotion_rad_s = sqrt( ...
        Earth.Mu_km3_s2 / OrbitCase.Radius_km^3);

    OrbitCase.Period_s = ...
        2*pi / OrbitCase.MeanMotion_rad_s;

    OrbitCase.Speed_km_s = sqrt( ...
        Earth.Mu_km3_s2 / OrbitCase.Radius_km);

    j2Coefficient = -1.5 * Earth.J2 * ...
        OrbitCase.MeanMotion_rad_s * ...
        (Earth.Radius_km/OrbitCase.Radius_km)^2;

    if strcmpi(OrbitCase.Type,"SSO")

        desiredRate_rad_s = ...
            2*pi/(tropicalYear_days*secondsPerDay);

        cosineInclination = ...
            desiredRate_rad_s/j2Coefficient;

        assert(abs(cosineInclination)<=1, ...
            'No first-order SSO solution exists at %.1f km.', ...
            altitude_km);

        OrbitCase.Inclination_deg = ...
            acosd(cosineInclination);
    end

    OrbitCase.RAANRate_rad_s = ...
        j2Coefficient*cosd(OrbitCase.Inclination_deg);

    OrbitCase.RAANRate_deg_day = ...
        rad2deg(OrbitCase.RAANRate_rad_s)*secondsPerDay;
end



function samples = buildDebrisOrientationSamples(spacecraftInclination_deg,Debris)
% Build deterministic circular-debris orbit orientation samples.
%
% The debris inclination distribution is sampled directly. Relative RAAN is
% uniformly distributed. For each pair of orbital planes, the mutual
% inclination Phi is calculated from:
%
% cos(Phi) = cos(i_s)cos(i_d) + sin(i_s)sin(i_d)cos(DeltaOmega)
%
% At an orbital-plane intersection, the debris velocity is approximated as
% lying in the spacecraft local T-N plane:
%
% v_d_hat = cos(Phi) T_hat +/- sin(Phi) N_hat
%
% The +/- sign represents the two possible node-crossing orientations.
% This avoids full debris propagation while preserving physically plausible
% circular-orbit velocity directions with zero radial velocity.

    inclination_deg = linspace(0,180,Debris.NumInclinationSamples);
    relativeRAAN_deg = linspace(0,360,Debris.NumRAANSamples+1);
    relativeRAAN_deg(end) = [];

    if strcmpi(Debris.InclinationModel,"Uniform")
        inclinationWeight = ones(size(inclination_deg));
    else
        inclinationWeight = Debris.InclinationFunction(inclination_deg);
    end

    assert(all(isfinite(inclinationWeight)) && ...
        all(inclinationWeight>=0), ...
        'Inclination weighting function must return finite nonnegative weights.');

    assert(sum(inclinationWeight)>0,'Inclination weights sum to zero.');
    inclinationWeight = inclinationWeight / sum(inclinationWeight);

    nSamples = Debris.NumInclinationSamples * ...
        Debris.NumRAANSamples * 2;

    samples.MutualInclination_deg = zeros(nSamples,1);
    samples.NodeSign = zeros(nSamples,1);
    samples.Weight = zeros(nSamples,1);

    index = 0;

    for iInc = 1:numel(inclination_deg)
        for iRAAN = 1:numel(relativeRAAN_deg)

            cosPhi = cosd(spacecraftInclination_deg)*cosd(inclination_deg(iInc)) + ...
                sind(spacecraftInclination_deg)*sind(inclination_deg(iInc))* ...
                cosd(relativeRAAN_deg(iRAAN));

            cosPhi = max(-1,min(1,cosPhi));
            phi_deg = acosd(cosPhi);

            for nodeSign = [-1 1]
                index = index + 1;

                samples.MutualInclination_deg(index) = phi_deg;
                samples.NodeSign(index) = nodeSign;

                samples.Weight(index) = ...
                    inclinationWeight(iInc) / ...
                    Debris.NumRAANSamples / 2;
            end
        end
    end

    samples.Weight = samples.Weight / sum(samples.Weight);
end



function Geometry = buildEncounterGeometry(Orbit,Shell,Sensor,Pointing,Debris,Exclusion,Numerics,Earth)
% Deterministic projected-area integration. No random objects or trajectories.
% For a uniform density and a constant velocity, integrating dA over lines
% intersecting the sensor gives its projected encounter area. Residence time
% is the observable chord length divided by relative speed on that line.
    samples = buildDebrisOrientationSamples(Orbit.Inclination_deg,Debris);
    % Combine numerically identical orientations without changing their total
    % probability. Original inclination/RAAN/node weights remain unchanged.
    [directions,~,group] = unique(round([samples.MutualInclination_deg ...
        samples.NodeSign],10),'rows');
    weights = accumarray(group,samples.Weight);
    edges = Shell.MinAltitude_km:Numerics.VelocityAltitudeStep_km:Shell.MaxAltitude_km;
    if edges(end)<Shell.MaxAltitude_km, edges(end+1)=Shell.MaxAltitude_km; end
    signs = 1;
    if Pointing.ChooseAntiSunNormal, signs = [-1 1]; end
    Geometry = cell(numel(signs),1);
    for g = 1:numel(signs)
        b = getBoresightRTN(Pointing,signs(g));
        % Square roll convention: one tangent axis is projected local T;
        % fallback to N if boresight is along T. Circular FOV ignores roll.
        a = [0 1 0]-b(2)*b;
        if norm(a)<1e-10, a = [0 0 1]-b(3)*b; end
        a = a/norm(a); c = cross(b,a);
        chunks = cell((numel(edges)-1)*numel(weights),1);
        index = 0;
        fprintf('Building FOV crossing geometry (normal sign %+d)...\n',signs(g));
        for ih = 1:numel(edges)-1
            h = (edges(ih)+edges(ih+1))/2;
            debrisSpeed = sqrt(Earth.Mu_km3_s2/(Earth.Radius_km+h));
            for j = 1:numel(weights)
                phi = directions(j,1); node = directions(j,2);
                velocity = debrisSpeed*[0 cosd(phi) node*sind(phi)]-[0 Orbit.Speed_km_s 0];
                speed = norm(velocity);
                index = index+1;
                if speed<1e-12 || weights(j)==0, continue; end
                v = velocity/speed;
                [p,dA] = projectedRays(b,v,Sensor,Numerics);
                [mid,L,ray] = observableChords(p,v,b,a,c,Orbit,Shell,Sensor,Exclusion,Earth);
                if isempty(L), continue; end
                altitude = vecnorm(mid+[Orbit.Radius_km 0 0],2,2)-Earth.Radius_km;
                select = altitude>=edges(ih) & altitude<edges(ih+1);
                if ih==numel(edges)-1, select = altitude>=edges(ih) & altitude<=edges(ih+1); end
                % Columns: RTN midpoint, residence seconds, rate coefficient
                % v*dA*orientation_weight [km^3/s], geocentric altitude.
                chunks{index} = [mid(select,:),L(select)/speed, ...
                    speed*dA(ray(select))*weights(j),altitude(select)];
            end
        end
        values = vertcat(chunks{:});
        if isempty(values), values = zeros(0,6); end
        Geometry{g}.Values = values;
        Geometry{g}.Boresight = b;
        Geometry{g}.NormalSign = signs(g);
        fprintf('  %d weighted observable chord samples retained.\n',size(values,1));
    end
end

function b = getBoresightRTN(Pointing,normalSign)
% RTN: R is radial out, T is along-track, N = R cross T.
    switch lower(Pointing.Mode)
        case "zenith", b = [1 0 0];
        case "alongtrack", b = [0 1 0];
        case "normal+", b = [0 0 1];
        case "normal-", b = [0 0 -1];
        case "customrtn", b = Pointing.CustomRTN;
        otherwise, error('Unknown pointing mode.');
    end
    validateattributes(b,{'numeric'},{'vector','numel',3,'real','finite'});
    b = b(:)'; assert(norm(b)>0,'Boresight cannot be zero.');
    if Pointing.ChooseAntiSunNormal, b(3)=normalSign*abs(b(3)); end
    b = b/norm(b);
end

function [p,dA] = projectedRays(b,v,Sensor,Numerics)
% Lines are parallel to relative velocity; their base points lie in the
% plane normal to v. A cubic spacing concentrates area samples near the apex
% so very short qualified ranges are not hidden by a uniform 500-km grid.
    R = Sensor.Range_km;
    beta = Sensor.FOV_deg/2;
    if strcmpi(Sensor.FOVShape,'Square'), beta=atand(sqrt(2)*tand(beta)); end
    e1 = b-dot(b,v)*v;
    if norm(e1)<1e-10
        e1 = cross(v,[1 0 0]);
        if norm(e1)<1e-10, e1=cross(v,[0 1 0]); end
    end
    e1 = e1/norm(e1); e2 = cross(v,e1);
    alpha = acosd(max(-1,min(1,dot(b,e1))));
    xmax = R*max(0,cosd(max(alpha-beta,0)));
    xmin = -R*max(0,cosd(max(180-alpha-beta,0)));
    n = Numerics.ProjectedAlongSamples;
    % Inverse signed cubic coordinate gives fine resolution around zero.
    u = linspace(sign(xmin)*abs(xmin)^(1/3),xmax^(1/3),n+1);
    xEdges = u.^3;
    x = (xEdges(1:end-1)+xEdges(2:end))/2;
    dx = diff(xEdges);
    sinFlowAngle = norm(cross(b,v));
    ymax = R*sind(beta)*ones(size(x));
    if sinFlowAngle>sind(beta)+1e-12
        slope = sind(beta)/sqrt(sinFlowAngle^2-sind(beta)^2);
        ymax = min(ymax,max(x,0)*slope);
    end
    ny = Numerics.ProjectedAcrossSamples;
    yFraction = ((1:ny)-0.5)/ny*2-1;
    [xx,yy] = ndgrid(x,yFraction);
    yy = yy.*ymax(:);
    p = xx(:)*e1+yy(:)*e2;
    area = repmat((dx.*2.*ymax/ny)',1,ny);
    dA = area(:);
end

function [mid,L,ray] = observableChords(p,v,b,a,c,Orbit,Shell,Sensor,Exclusion,Earth)
% Find all line intersections with physical boundaries; test each interval.
% Merge neighboring valid intervals so irrelevant mathematical roots do not
% create multiple detections. Shell re-entry across a real excluded gap is a
% separate observable passage. Earth-limb masking also prevents occultation.
    n = size(p,1); R = Sensor.Range_km;
    bound = [-R*ones(n,1),R*ones(n,1),sphereRoots(p,v,R)];
    if Sensor.MinRange_km>0
        bound = [bound sphereRoots(p,v,Sensor.MinRange_km)];
    end
    if strcmpi(Sensor.FOVShape,'Circular')
        bound = [bound coneRoots(p,v,b,Sensor.FOV_deg/2)];
    else
        k = tand(Sensor.FOV_deg/2);
        for normal = [a-k*b;-a-k*b;c-k*b;-c-k*b]'
            denom = v*normal;
            if abs(denom)>1e-14, bound=[bound -(p*normal)/denom]; end %#ok<AGROW>
        end
    end
    if abs(dot(v,b))>1e-14, bound=[bound -(p*b')/dot(v,b)]; end
    q = p+[Orbit.Radius_km 0 0];
    bound = [bound sphereRoots(q,v,Earth.Radius_km+Shell.MinAltitude_km) ...
        sphereRoots(q,v,Earth.Radius_km+Shell.MaxAltitude_km)];
    limb = asind(Earth.Radius_km/Orbit.Radius_km)+Exclusion.EarthLimbClearance_deg;
    assert(limb<90,'Earth margin is too large for this orbit.');
    bound = [bound coneRoots(p,v,[-1 0 0],limb)];
    bound(~isfinite(bound)) = Inf;
    bound = sort(bound,2);
    duplicate = [false(n,1),diff(bound,1,2)<1e-9];
    bound(duplicate) = Inf;
    bound = sort(bound,2);
    t0 = bound(:,1:end-1); t1 = bound(:,2:end);
    t = (t0+t1)/2;
    valid = isfinite(t) & t1-t0>1e-10 & t>=-R & t<=R;
    for j = 1:size(t,2)
        point = p+t(:,j)*v;
        range = vecnorm(point,2,2);
        height = vecnorm(point+[Orbit.Radius_km 0 0],2,2)-Earth.Radius_km;
        along = point*b';
        if strcmpi(Sensor.FOVShape,'Circular')
            inFOV = along>=range*cosd(Sensor.FOV_deg/2);
        else
            inFOV = along>=0 & abs(point*a')<=along*tand(Sensor.FOV_deg/2) & ...
                abs(point*c')<=along*tand(Sensor.FOV_deg/2);
        end
        clearEarth = -point(:,1)<range*cosd(limb);
        valid(:,j) = valid(:,j) & range>0 & range>=Sensor.MinRange_km & ...
            range<=R & height>=Shell.MinAltitude_km & height<=Shell.MaxAltitude_km & ...
            inFOV & clearEarth;
    end
    starts = [valid(:,1),valid(:,2:end)&~valid(:,1:end-1)];
    ends = [valid(:,1:end-1)&~valid(:,2:end),valid(:,end)];
    % Transpose before find so starts and ends are paired in ray-major order.
    [startCol,ray] = find(starts');
    [endCol,endRay] = find(ends');
    assert(isequal(ray,endRay),'Internal chord pairing error.');
    first = t0(sub2ind(size(t0),ray,startCol));
    last = t1(sub2ind(size(t1),ray,endCol));
    L = last-first;
    mid = p(ray,:)+(first+last)/2*v;
end

function roots = sphereRoots(p,v,radius)
    B = p*v';
    discriminant = B.^2-sum(p.^2,2)+radius^2;
    s = sqrt(max(discriminant,0));
    roots = [-B-s,-B+s];
    roots(discriminant<0,:) = NaN;
end

function roots = coneRoots(p,v,axis,halfAngle)
% The quadratic includes both cone nappes. The interval validity test selects
% the correct forward cone, so the opposite nappe cannot admit detections.
    cs = cosd(halfAngle)^2;
    pa = p*axis'; va = dot(v,axis);
    A = va^2-cs;
    B = 2*(pa*va-cs*(p*v'));
    C = pa.^2-cs*sum(p.^2,2);
    roots = nan(size(p,1),2);
    if abs(A)<1e-13
        good = abs(B)>1e-13;
        roots(good,1) = -C(good)./B(good);
    else
        disc = B.^2-4*A*C;
        good = disc>=0;
        roots(good,:) = [(-B(good)-sqrt(disc(good)))/(2*A), ...
            (-B(good)+sqrt(disc(good)))/(2*A)];
    end
end

function [rates,Diagnostics] = evaluatePlanes(offsets,Geometry,Orbit,Mission,Sensor,Response, ...
    Environment,Debris,Pointing,Exclusion,Numerics,Earth)
% Solar direction uses a low-order mean-Sun/ecliptic model (no toolbox).
% SSO node tracks mean Sun RA at the requested LTAN. For non-SSO Circular
% mode, initialize the same LTAN and use first-order J2 node precession.
% Neither construction is a precision orbital ephemeris.
    nBins = size(Debris.SizeBins_cm,1);
    rates = zeros(numel(offsets),nBins);
    Diagnostics.MinPhase_deg = NaN; Diagnostics.MaxPhase_deg = NaN;
    phaseMin = Inf; phaseMax = -Inf;
    nDates = Numerics.MissionDateSamples;
    nOrbit = Numerics.OrbitalPositionSamples;
    beta = Sensor.FOV_deg/2;
    if strcmpi(Sensor.FOVShape,'Square'), beta=atand(sqrt(2)*tand(beta)); end
    [~,sunRA0] = approximateSun(Mission.StartDate);
    for plane = 1:numel(offsets)
        for d = 1:nDates
            elapsed = (d-0.5)/nDates*Mission.Duration_days;
            date = Mission.StartDate+days(elapsed);
            [sunECI,sunRA] = approximateSun(date);
            ltan = Orbit.NominalLTAN_hr+offsets(plane)/60;
            if strcmpi(Orbit.Type,'SSO')
                raan = sunRA+15*(ltan-12);
            else
                raan = sunRA0+15*(ltan-12)+Orbit.RAANRate_deg_day*elapsed;
            end
            for u = ((1:nOrbit)-0.5)/nOrbit*360
                Q = rtnToECI(raan,Orbit.Inclination_deg,u);
                sunRTN = sunECI*Q;
                g = 1;
                if Pointing.ChooseAntiSunNormal
                    signN = 1; if sunRTN(3)>0, signN=-1; end
                    g = find(cellfun(@(x) x.NormalSign==signN,Geometry),1);
                end
                G = Geometry{g};
                sunSeparation = acosd(max(-1,min(1,dot(G.Boresight,sunRTN))));
                % Conservative solar angular radius at approximately perihelion.
                if sunSeparation<=beta+Exclusion.SunClearance_deg+0.272, continue; end
                V = G.Values;
                if isempty(V), continue; end
                mid = V(:,1:3);
                range = vecnorm(mid,2,2);
                phase = [];
                if Sensor.UsePhaseAngleCSV || Sensor.ReportPhaseDiagnostic
                    % Angle at debris: direction toward sensor is -mid; Sun
                    % is at infinity here (parallax negligible at this fidelity).
                    phase = acosd(max(-1,min(1,(-mid*sunRTN')./range)));
                    phaseMin = min(phaseMin,min(phase)); phaseMax = max(phaseMax,max(phase));
                end
                for bin = 1:nBins
                    diameter = Debris.RepresentativeDiameter_cm(bin)*ones(size(range));
                    density = queryDensity(Environment,bin,V(:,6));
                    if Sensor.UsePhaseAngleCSV
                        needed = Response.Time(range,diameter,phase);
                    else
                        needed = Response.Time(range,diameter);
                    end
                    if any(~isfinite(needed) & density>0)
                        error(['Sensor CSV coverage missing at a sampled range/diameter/phase. ' ...
                            'No extrapolation was performed. Supply a wider response table.']);
                    end
                    qualified = V(:,4)>=needed & isfinite(needed);
                    rates(plane,bin) = rates(plane,bin)+ ...
                        sum(density(qualified).*V(qualified,5))/(nDates*nOrbit);
                end
            end
        end
    end
    if isfinite(phaseMin)
        Diagnostics.MinPhase_deg = phaseMin; Diagnostics.MaxPhase_deg = phaseMax;
    end
    Diagnostics.DateSamples = nDates; Diagnostics.OrbitSamples = nOrbit;
    Diagnostics.Note = 'Phases describe sampled observable chords, not only qualified detections.';
end

function Q = rtnToECI(raan,inclination,u)
% Columns are R,T,N expressed in ECI; row RTN vectors transform with Q'.
    p = [cosd(raan);sind(raan);0];
    q = [-sind(raan)*cosd(inclination);cosd(raan)*cosd(inclination);sind(inclination)];
    R = p*cosd(u)+q*sind(u);
    T = -p*sind(u)+q*cosd(u);
    Q = [R T cross(R,T)];
end

function [sun,meanRA] = approximateSun(date)
% Low-order geometric solar direction relative to J2000; mean solar longitude
% advances uniformly. Adequate for architecture sampling, not precise margins.
    n = days(date-datetime(2000,1,1,12,0,0));
    meanLongitude = mod(280.460+0.9856474*n,360);
    anomaly = mod(357.528+0.9856003*n,360);
    lambda = meanLongitude+1.915*sind(anomaly)+0.020*sind(2*anomaly);
    epsilon = 23.439-0.0000004*n;
    sun = [cosd(lambda),cosd(epsilon)*sind(lambda),sind(epsilon)*sind(lambda)];
    sun = sun/norm(sun);
    % Mean-Sun RA convention for local MEAN solar time. The orbit node uses
    % this uniformly advancing angle; exclusion uses the geometric Sun above.
    meanRA = meanLongitude;
end

function offsets = distributePlanesInLTAN(numPlanes,totalSpread_min)
    if numPlanes==1, offsets=0; else, offsets=linspace(-totalSpread_min/2,totalSpread_min/2,numPlanes); end
end

function plotTradeResults(Architecture,LTAN,Altitude,runLabel)
    if ~isempty(Architecture)
        figure('Color','w','Name','Architecture trade');
        p = unique(Architecture.Planes); n = unique(Architecture.SpacecraftPerPlane);
        Z = reshape(Architecture.TotalPerDay,numel(n),numel(p));
        imagesc(p,n,Z); axis xy; colorbar;
        xlabel('Orbital planes'); ylabel('Spacecraft per plane');
        title({'Constellation qualified passages / day',char(runLabel)});
    end
    if ~isempty(LTAN)
        figure('Color','w','Name','LTAN trade');
        plot(LTAN.LTANSpread_min,[LTAN.SciencePerDay LTAN.OutsidePerDay],'LineWidth',1.5);
        xlabel('Total LTAN spread (min)'); ylabel('Constellation passages / day');
        legend('Science sizes','Outside science sizes'); grid on;
        title(char(runLabel));
    end
    if ~isempty(Altitude)
        figure('Color','w','Name','Altitude trade');
        plot(Altitude.Altitude_km,[Altitude.SciencePerDay Altitude.OutsidePerDay],'LineWidth',1.5);
        xlabel('Spacecraft altitude (km)'); ylabel('Constellation passages / day');
        legend('Science sizes','Outside science sizes'); grid on;
        title(char(runLabel));
    end
end

function runGeometrySelfChecks
% Small analytic cases catch accidental changes to chord/area conventions.
% These do not establish convergence or scientific validity of the mission run.
    Sensor.Range_km = 10; Sensor.MinRange_km = 0; Sensor.FOV_deg = 20;
    Sensor.FOVShape = "Circular";
    Numerics.ProjectedAlongSamples = 160; Numerics.ProjectedAcrossSamples = 80;
    Orbit.Radius_km = 7000;
    Shell.MinAltitude_km = 1; Shell.MaxAltitude_km = 20000;
    Exclusion.EarthLimbClearance_deg = 0.05; Earth.Radius_km = 6371;
    b = [1 0 0]; a = [0 1 0]; c = [0 0 1];
    for v = [0 1 0;1 0 0]'
        velocityDirection = v';
        [p,dA] = projectedRays(b,velocityDirection,Sensor,Numerics);
        [~,L,ray] = observableChords(p,velocityDirection,b,a,c,Orbit,Shell,Sensor,Exclusion,Earth);
        if abs(dot(b,velocityDirection))<0.5
            exact = Sensor.Range_km^2*deg2rad(Sensor.FOV_deg/2);
        else
            exact = pi*(Sensor.Range_km*sind(Sensor.FOV_deg/2))^2;
        end
        assert(all(L>0) && abs(sum(dA(ray))/exact-1)<0.012, ...
            'Projected-area self-check failed.');
    end
    b = [-1 0 0]; c = cross(b,a);
    [p,~] = projectedRays(b,[0 1 0],Sensor,Numerics);
    [~,L,~] = observableChords(p,[0 1 0],b,a,c,Orbit,Shell,Sensor,Exclusion,Earth);
    assert(isempty(L),'Earth-background exclusion self-check failed.');
    fprintf('Geometry self-checks passed.\n');
end
