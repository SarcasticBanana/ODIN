function [rate_per_s,details] = debrisDetectionRate(time_s,r_km,v_km_s,h_km,cfg,sunPosition_km,Detection)
% DEBRISDETECTIONRATE Qualified passage-rate history for one spacecraft.
% Updated 2026-10-29: direct CDH Cartesian time histories, one-file interface.
%
% [rates,details] = debrisDetectionRate(t,r,v,h,cfg,rSun,Detection)
% Detection is optional; grouped defaults follow below.
%
% INPUTS (all rows describe the SAME times and fixed inertial frame)
%   time_s         N-vector of finite, strictly increasing sample times [s].
%                  One sample is also supported. No resampling is performed.
%   r_km           N-by-3 Earth-centered spacecraft POSITION [km].
%   v_km_s         N-by-3 inertial spacecraft VELOCITY [km/s], not ECEF or
%                  velocity relative to the rotating atmosphere.
%   h_km           N-vector CDH altitude [km], NOT angular momentum or radius.
%                  Preserved for comparison with |r|-cfg.Re. Geometry uses r.
%                  Ellipsoidal CDH drag height can differ from spherical
%                  MASTER altitude; differences >1 metre generate one warning.
%   cfg            CDH constants: mu [km^3/s^2], Re [km], Rsun [km], epoch.
%                  epoch identifies the caller's frame; it does not generate
%                  Sun positions or propagate states at the supplied times.
%   sunPosition_km N-by-3 simultaneous Earth-to-Sun POSITION [km], same frame
%                  as r and v. No unit directions, implicit Sun broadcasting,
%                  or internal Sun ephemeris. Supply the CDH Sun history.
%   Detection      Optional scalar configuration; same settings as before.
%
% Supply N-by-3 arrays (transpose CDH 3-by-N outputs explicitly). For a single
% time, Cartesian vectors can be rows or columns. Scalars r and v containing
% only radius and speed are insufficient to establish pointing and exclusion.
% This interface replaces the earlier orbital-elements calling signature.
%
% OUTPUTS
%   rate_per_s     N-by-B, expected qualified OBJECT PASSAGES per second.
%                  Rows align exactly with time_s; B is the number of bins.
%                  Default bins: [1,1.5], [1.5,2.5], [2.5,3], [3,5] cm.
%   details        Time histories of state, pointing, exclusion, qualification
%                  diagnostics and total/science/outside-science rates; shared
%                  MASTER/sensor metadata stored once. SizeBins_cm defines
%                  output columns. GeometryBuildCount reports cache misses.
%
% EXAMPLE (variables supplied by the CDH orbit and Sun models)
%   Detection.Debris.MASTERFolder = fullfile(pwd,'MASTER_Data');
%   Detection.Sensor.CSVFile = fullfile(pwd,'Sensor_Residence_Time_Trade.csv');
%   [rate_per_s,info] = debrisDetectionRate(t,r,v,h,cfg,rSun,Detection);
%   % Precompute before the CDH ODE. Per-bin lookup, within the sampled span:
%   % lambda = interp1(t(:),rate_per_s,tQuery,'linear');
%   % Do not interpolate blindly across exclusion/operating transitions:
%   % resolve transition times and enforce the observing mask in the caller.
%
% PHYSICS AND CDH CONVENTIONS
%   Retains MASTER altitude-diameter density, circular debris inclination/RAAN
%   sampling, sensor residence-time CSV, pointing, Earth and Sun exclusion.
%   Spacecraft motion is not assumed circular: full radial/transverse velocity
%   comes from r,v. MASTER queries use spherical altitude, not CDH drag height.
%   RTN is R=r/|r|, N=(r cross v)/|r cross v|, T=N cross R. Inclination comes
%   from N relative to the fixed frame's equatorial plane.
%   CustomRTN is fixed in that local frame across the input history, with the
%   existing optional anti-Sun normal-sign choice. No new pointing schedule.
%   Illumination remains in the CSV; fixed-phase limitations remain. No added
%   spacecraft eclipse multiplier, velocity threshold, or brightness rescaling.
%   Straight chords, locally fixed attitude/velocity and midpoint response/
%   density remain approximations. No individual debris propagation.
%
% PERFORMANCE / SCOPE
%   MASTER and sensor files are loaded/checked once per batch (and cached
%   between calls). Exact repeated LAST geometry is reused; distinct states
%   still require chord integration. No approximate state rounding or physics
%   change is introduced. A history interface alone does not guarantee speedup.
%   No plots, trades, comms downtime, synthetic population, constellation,
%   automatic saving, count integration or ODE calls. Caller owns these.
%   After replacing files while preserving timestamps: clear debrisDetectionRate
%   Use a configurable CDH sampling cadence and check rate/count convergence.
%
%% 1. DETECTION INPUTS AND DEFAULTS

    if nargin<7 || isempty(Detection), Detection=struct; end
    modelFolder = fileparts(mfilename('fullpath'));
    Defaults.Shell.MinAltitude_km = 600;
    Defaults.Shell.MaxAltitude_km = 1000;

    Defaults.Sensor.CSVFile = fullfile(modelFolder,'Sensor_Residence_Time_Trade.csv');
    Defaults.Sensor.FOVShape = "Square"; % "Circular" or "Square"
    Defaults.Sensor.UseCSVFOV = true;
    Defaults.Sensor.FOV_deg = 1.4;         % FULL angle if CSV FOV is disabled
    Defaults.Sensor.ExpectedAlbedo = 0.20;% Metadata check, not brightness scaling
    Defaults.Sensor.UsePhaseAngleCSV = false;
    Defaults.Sensor.FixedPhaseAngle_deg = 30;
    Defaults.Sensor.ReportPhaseDiagnostic = false;

    Defaults.Pointing.Mode = "CustomRTN"; % Zenith, AlongTrack, Velocity, Normal+/-, CustomRTN
    Defaults.Pointing.CustomRTN = [sind(66) 0 cosd(66)];
    Defaults.Pointing.ChooseAntiSunNormal = true;
    % R is outward; N = r cross v; T = N cross R. For eccentric orbits T is
    % transverse, not exactly parallel to velocity. "AlongTrack" retains +T;
    % "Velocity" follows the actual spacecraft velocity including radial motion.
    % Square FOV roll uses local T projected normal to boresight (N fallback).

    Defaults.Debris.MASTERFolder = fullfile(modelFolder,'MASTER_Data');
    Defaults.Debris.MASTERFileName = "ODIN_2024_D010_070_t.__1";
    Defaults.Debris.SizeBins_cm = [1 1.5;1.5 2.5;2.5 3;3 5];
    Defaults.Debris.RepresentativeDiameter_cm = [1.25 2.5 3 5];
    Defaults.Debris.ScienceSizeRange_cm = [1 3];
    Defaults.Debris.MaxUsableDiameter_cm = 5;
    Defaults.Debris.ExpectedReferenceEpoch = "";
    Defaults.Debris.InclinationModel = "Uniform";
    Defaults.Debris.InclinationFunction = @(i_deg) ones(size(i_deg));
    Defaults.Debris.NumInclinationSamples = 100;
    Defaults.Debris.NumRAANSamples = 8;
    % All bin edges must align with MASTER diameter cells. Only configured
    % bins are counted, not the extra upper-diameter buffer in the export.

    Defaults.Exclusion.SunClearance_deg = 5;
    Defaults.Exclusion.EarthLimbClearance_deg = 0.05;
    % Sun clearance is extra angle beyond the apparent solar disk and FOV.
    % Square FOV uses a conservative circumscribed cone for Sun exclusion.

    Defaults.Numerics.VelocityAltitudeStep_km = 50;
    Defaults.Numerics.ProjectedAlongSamples = 40;
    Defaults.Numerics.ProjectedAcrossSamples = 12;

    Detection = mergeOptions(Defaults,Detection,'Detection');
    Shell = Detection.Shell; Sensor = Detection.Sensor;
    Debris = Detection.Debris; Pointing = Detection.Pointing;
    Exclusion = Detection.Exclusion; Numerics = Detection.Numerics;
    validateDetectionInputs(Detection);

%% 2. VALIDATE CDH TIME HISTORY

    assert(nargin>=6,'Provide time_s, r_km, v_km_s, h_km, cfg, and sunPosition_km.');
    validateattributes(time_s,{'numeric'},{'vector','nonempty','real','finite'});
    time_s = time_s(:);
    nTimes = numel(time_s);
    assert(all(diff(time_s)>0),'time_s must be strictly increasing, in seconds.');
    r_km = validateHistoryVectors(r_km,nTimes,'r_km');
    v_km_s = validateHistoryVectors(v_km_s,nTimes,'v_km_s');
    sunPosition_km = validateHistoryVectors(sunPosition_km,nTimes,'sunPosition_km');
    validateattributes(h_km,{'numeric'},{'vector','numel',nTimes,'real','finite'});
    h_km = h_km(:);
    assert(isstruct(cfg) && isscalar(cfg),'cfg must be a scalar CDH configuration struct.');
    for field = {'mu','Re','Rsun'}
        name = field{1};
        assert(isfield(cfg,name),'Missing CDH constant cfg.%s.',name);
        validateattributes(cfg.(name),{'numeric'},{'scalar','real','finite','positive'});
    end
    assert(isfield(cfg,'epoch') && ~isempty(cfg.epoch), ...
        'cfg.epoch is needed to identify the caller''s fixed reference frame.');
    Earth.Mu_km3_s2 = cfg.mu;
    Earth.Radius_km = cfg.Re;

    % Position is authoritative for geometry and MASTER's spherical altitude.
    % h_km is preserved as CDH altitude metadata: CDH may use ellipsoidal
    % drag height, which must not silently replace spherical MASTER altitude.
    sphericalAltitude_km = vecnorm(r_km,2,2)-cfg.Re;
    assert(all(sphericalAltitude_km>0),'A spacecraft state is inside the Earth sphere.');
    altitudeDifference_km = h_km-sphericalAltitude_km;
    if any(abs(altitudeDifference_km)>0.001)
        warning('Detection:AltitudeConvention', ...
            ['CDH h differs from |r|-cfg.Re by up to %.6g km. Position r controls ' ...
             'geometry; h is retained in details. Check altitude convention and units.'], ...
            max(abs(altitudeDifference_km)));
    end

%% 3. LOAD AND CACHE MASTER DENSITY / SENSOR CSV

    persistent cachedData cachedGeometry
    masterPath = locateMASTERFile(Debris,Shell);
    Debris.MASTERFolder = fileparts(masterPath);
    [~,stem,extension] = fileparts(masterPath);
    Debris.MASTERFileName = string([stem extension]);
    key = {fileStamp(masterPath),fileStamp(Sensor.CSVFile),Shell, ...
        Debris.SizeBins_cm,Debris.MaxUsableDiameter_cm, ...
        Debris.ExpectedReferenceEpoch,Sensor,Earth.Radius_km};
    if isempty(cachedData) || ~isequaln(key,cachedData.Key)
        [Response,Sensor] = readSensorResponse(Sensor,Debris);
        Environment = readMASTEREnvironment(Debris,Shell,Earth);
        cachedData = struct('Key',{key},'Response',Response, ...
            'Sensor',Sensor,'Environment',Environment);
    else
        Response = cachedData.Response;
        Sensor = cachedData.Sensor;
        Environment = cachedData.Environment;
        assert(all(Debris.RepresentativeDiameter_cm>=min(Response.Diameter_cm) & ...
            Debris.RepresentativeDiameter_cm<=max(Response.Diameter_cm)), ...
            'Representative diameters exceed sensor CSV coverage.');
    end

%% 4. EVALUATE EACH SUPPLIED TIME (NO ORBIT PROPAGATION)

    historyRates = zeros(nTimes,size(Debris.SizeBins_cm,1));
    historyStatus = strings(nTimes,1);
    historyBoresight = zeros(nTimes,3);
    historyVelocityRTN = zeros(nTimes,3);
    historySunSeparation = zeros(nTimes,1);
    historySunRadius = zeros(nTimes,1);
    historySunExcluded = false(nTimes,1);
    historyPhase = nan(nTimes,2);
    historyChordCount = zeros(nTimes,1);
    historyQualified = zeros(size(historyRates));
    historyInclination = zeros(nTimes,1);
    geometryBuilds = 0;
    % figure;
    for it = 1:nTimes
        % disp(['looking at step ', num2str(it)])
        [Orbit,Q,rSpacecraft_km,vSpacecraft_km_s] = ...
            spacecraftFromCartesian(r_km(it,:),v_km_s(it,:),Earth);
        sunAtTime_km = sunPosition_km(it,:);
        validateattributes(sunAtTime_km,{'numeric'},{'vector','numel',3,'real','finite'});
        sunAtTime_km = sunAtTime_km(:)';
        sunFromSpacecraft_km = sunAtTime_km-rSpacecraft_km;
        sunDistance_km = norm(sunFromSpacecraft_km);
        assert(sunDistance_km>cfg.Rsun && norm(sunAtTime_km)>cfg.Rsun, ...
            'Sun input must be an Earth-centered position in km, not a unit direction.');
        sunRTN_km = sunFromSpacecraft_km*Q;
        sunRTN = sunRTN_km/sunDistance_km;
        Pointing.CurrentVelocityRTN_km_s = Orbit.VelocityRTN_km_s;
        normalSign = 1;
        if sunRTN(3)>0, normalSign=-1; end
        boresight = getBoresightRTN(Pointing,normalSign);
        sunSeparation_deg = acosd(max(-1,min(1,dot(boresight,sunRTN))));
        sunRadius_deg = asind(cfg.Rsun/sunDistance_km);

        % Existing instantaneous qualification, unchanged for each state.
        nBins = size(Debris.SizeBins_cm,1);
        rate_per_s = zeros(1,nBins);
        science = Debris.SizeBins_cm(:,1)>=Debris.ScienceSizeRange_cm(1) & ...
            Debris.SizeBins_cm(:,2)<=Debris.ScienceSizeRange_cm(2);
        beta_deg = Sensor.FOV_deg/2;
        if strcmpi(Sensor.FOVShape,'Square'), beta_deg=atand(sqrt(2)*tand(beta_deg)); end
        sunBlocked = sunSeparation_deg<=beta_deg+Exclusion.SunClearance_deg+sunRadius_deg;

        details.Status = "NoQualifiedPassages";
        details.SizeBins_cm = Debris.SizeBins_cm;
        details.RepresentativeDiameter_cm = Debris.RepresentativeDiameter_cm;
        details.ScienceBin = science';
        details.Inclination_deg = Orbit.Inclination_deg;
        details.Position_km = rSpacecraft_km;
        details.Velocity_km_s = vSpacecraft_km_s;
        details.VelocityRTN_km_s = Orbit.VelocityRTN_km_s;
        details.SphericalAltitude_km = Orbit.Radius_km-cfg.Re;
        details.FrameEpoch = cfg.epoch;
        details.FrameConvention = "Caller mean equator/equinox of epoch; Sun must match";
        details.BoresightRTN = boresight;
        details.SunSeparation_deg = sunSeparation_deg;
        details.SunAngularRadius_deg = sunRadius_deg;
        details.SunExcluded = sunBlocked;
        details.PhaseRange_deg = [NaN NaN];
        details.PhaseInterpolationEnabled = Sensor.UsePhaseAngleCSV;
        details.FixedResponsePhase_deg = Sensor.FixedPhaseAngle_deg;
        details.FixedPhaseApproximation = ~Sensor.UsePhaseAngleCSV;
        details.SensorCSV = Sensor.CSVFile;
        details.SensorFOV_deg = Sensor.FOV_deg;
        details.MASTERFile = masterPath;
        details.MASTERReferenceEpoch = Environment.ReferenceEpoch;
        details.MASTERSimulationInterval = Environment.SimulationInterval;
        details.PopulationSummary = Environment.Summary;
        details.ObservableChordSamples = 0;
        details.QualifiedChordSamples = zeros(1,nBins);
        details.Numerics = Numerics;
        details.EarthExclusionModel = "Conservative equatorial-radius sphere plus limb margin";
        details.RateMeaning = "Qualified object passages per second for one spacecraft";

        if sunBlocked
            details.Status = "SunExcluded";
        else
            % Geometry is built only for the boresight used at this instant.
            effectivePointing = Pointing;
            effectivePointing.Mode = "CustomRTN";
            effectivePointing.CustomRTN = boresight;
            effectivePointing.ChooseAntiSunNormal = false;
            inclinationKey = {Debris.InclinationModel};
            if strcmpi(Debris.InclinationModel,'Function')
                inclinationKey{2} = Debris.InclinationFunction;
            end
            geometryKey = {Orbit.Radius_km,Orbit.Inclination_deg,Orbit.VelocityRTN_km_s, ...
                boresight,Shell,Sensor.Range_km,Sensor.MinRange_km,Sensor.FOV_deg, ...
                Sensor.FOVShape,inclinationKey, ...
                Debris.NumInclinationSamples,Debris.NumRAANSamples,Exclusion,Numerics,Earth};
            if isempty(cachedGeometry) || ~isequaln(geometryKey,cachedGeometry.Key)
                geometryBuilds = geometryBuilds+1;
                G = buildEncounterGeometry(Orbit,Shell,Sensor,effectivePointing,Debris, ...
                    Exclusion,Numerics,Earth);
                V = G{1}.Values;
                cachedGeometry = struct('Key',{geometryKey},'Values',V);
            else
                V = cachedGeometry.Values;
            end
            details.ObservableChordSamples = size(V,1);
            if isempty(V)
                details.Status = "NoObservableShellChords";
            else
                mid = V(:,1:3);
                range_km = vecnorm(mid,2,2);
                phase_deg = [];
                if Sensor.UsePhaseAngleCSV || Sensor.ReportPhaseDiagnostic
                    % At debris: target-to-sensor=-mid, target-to-Sun=sunRTN_km-mid.
                    targetToSun = sunRTN_km-mid;
                    phaseCos = sum((-mid).*targetToSun,2)./(range_km.*vecnorm(targetToSun,2,2));
                    phase_deg = acosd(max(-1,min(1,phaseCos)));
                    details.PhaseRange_deg = [min(phase_deg),max(phase_deg)];
                end
                for bin = 1:nBins
                    diameter = Debris.RepresentativeDiameter_cm(bin)*ones(size(range_km));
                    density = queryDensity(Environment,bin,V(:,6));
                    if Sensor.UsePhaseAngleCSV
                        needed_s = Response.Time(range_km,diameter,phase_deg);
                    else
                        needed_s = Response.Time(range_km,diameter);
                    end
                    assert(~any(~isfinite(needed_s) & density>0), ...
                        'Sensor CSV coverage is missing for a sampled range/diameter/phase; no extrapolation is allowed.');
                    qualified = isfinite(needed_s) & V(:,4)>=needed_s;
                    rate_per_s(bin) = sum(density(qualified).*V(qualified,5));
                    details.QualifiedChordSamples(bin) = nnz(qualified);
                end
                if any(rate_per_s>0), details.Status="OK"; end
            end
        end
        details.Rate_per_s = rate_per_s;
        details.ScienceRate_per_s = sum(rate_per_s(science));
        details.OutsideScienceRate_per_s = sum(rate_per_s(~science));
        details.TotalRate_per_s = sum(rate_per_s);
        historyRates(it,:) = rate_per_s;
        % plot(time_s(1:it), sum(historyRates(1:it,:),2),'o-');
        % drawnow;

        historyStatus(it) = details.Status;
        historyBoresight(it,:) = boresight;
        historyVelocityRTN(it,:) = Orbit.VelocityRTN_km_s;
        historySunSeparation(it) = sunSeparation_deg;
        historySunRadius(it) = sunRadius_deg;
        historySunExcluded(it) = sunBlocked;
        historyPhase(it,:) = details.PhaseRange_deg;
        historyChordCount(it) = details.ObservableChordSamples;
        historyQualified(it,:) = details.QualifiedChordSamples;
        historyInclination(it) = Orbit.Inclination_deg;
    end

%% 5. RETURN ALIGNED TIME HISTORIES AND SHARED INPUT METADATA

    % Keep shared sensor/environment metadata once, not once per sample.
    rate_per_s = historyRates;
    details.Time_s = time_s;
    details.Rate_per_s = historyRates;
    details.Status = historyStatus;
    details.Position_km = r_km;
    details.Velocity_km_s = v_km_s;
    details.InputAltitude_km = h_km;
    details.SphericalAltitude_km = sphericalAltitude_km;
    details.AltitudeDifference_km = altitudeDifference_km;
    details.AltitudePolicy = "Position controls geometry; supplied h is metadata";
    details.SunPosition_km = sunPosition_km;
    details.Inclination_deg = historyInclination;
    details.VelocityRTN_km_s = historyVelocityRTN;
    details.BoresightRTN = historyBoresight;
    details.SunSeparation_deg = historySunSeparation;
    details.SunAngularRadius_deg = historySunRadius;
    details.SunExcluded = historySunExcluded;
    details.PhaseRange_deg = historyPhase;
    details.ObservableChordSamples = historyChordCount;
    details.QualifiedChordSamples = historyQualified;
    details.ScienceRate_per_s = sum(historyRates(:,science),2);
    details.OutsideScienceRate_per_s = sum(historyRates(:,~science),2);
    details.TotalRate_per_s = sum(historyRates,2);
    details.GeometryBuildCount = geometryBuilds;
    details.RateMeaning = "One row per supplied time; qualified passages/s per size bin";
end

%% 6. LOCAL INPUT AND SPACECRAFT HELPERS

function output = mergeOptions(defaults,provided,prefix)
% Reject misspelled options rather than silently ignoring caller settings.
    assert(isstruct(provided) && isscalar(provided),'%s must be a scalar struct.',prefix);
    output = defaults;
    names = fieldnames(provided);
    for k = 1:numel(names)
        name = names{k};
        assert(isfield(defaults,name),'Unknown option %s.%s.',prefix,name);
        if isstruct(defaults.(name))
            output.(name) = mergeOptions(defaults.(name),provided.(name),[prefix '.' name]);
        else
            output.(name) = provided.(name);
        end
    end
end

function validateDetectionInputs(D)
    assert(D.Shell.MinAltitude_km>0 && D.Shell.MaxAltitude_km>D.Shell.MinAltitude_km);
    b = D.Debris.SizeBins_cm;
    validateattributes(b,{'numeric'},{'2d','ncols',2,'nonempty','real','finite','positive'});
    assert(all(b(:,2)>b(:,1)) && all(b(2:end,1)>=b(1:end-1,2)), ...
        'Size bins must be sorted, positive-width and disjoint.');
    reps = D.Debris.RepresentativeDiameter_cm;
    validateattributes(reps,{'numeric'},{'vector','real','finite','positive'});
    assert(numel(reps)==size(b,1) && all(reps(:)>=b(:,1) & reps(:)<=b(:,2)), ...
        'Supply one representative diameter within each configured bin.');
    assert(all(b(:,2)<=D.Debris.MaxUsableDiameter_cm),'Bin exceeds MaxUsableDiameter_cm.');
    science = D.Debris.ScienceSizeRange_cm;
    validateattributes(science,{'numeric'},{'vector','numel',2,'real','finite','positive'});
    assert(science(2)>science(1));
    inside = b(:,1)>=science(1) & b(:,2)<=science(2);
    straddle = b(:,1)<science(2) & b(:,2)>science(1) & ~inside;
    assert(~any(straddle),'Split size bins at science-range boundaries.');
    assert(any(strcmpi(D.Debris.InclinationModel,{'Uniform','Function'})), ...
        'InclinationModel must be Uniform or Function.');
    assert(isa(D.Debris.InclinationFunction,'function_handle'));
    counts = [D.Debris.NumInclinationSamples D.Debris.NumRAANSamples ...
        D.Numerics.ProjectedAlongSamples D.Numerics.ProjectedAcrossSamples];
    assert(all(isfinite(counts) & counts>=2 & mod(counts,1)==0));
    validateattributes(D.Numerics.VelocityAltitudeStep_km,{'numeric'},{'scalar','finite','positive'});
    margins = [D.Exclusion.SunClearance_deg D.Exclusion.EarthLimbClearance_deg];
    assert(all(isfinite(margins) & margins>=0),'Exclusion margins must be finite and nonnegative.');
    assert(isscalar(D.Sensor.UseCSVFOV) && isscalar(D.Sensor.UsePhaseAngleCSV) && ...
        isscalar(D.Sensor.ReportPhaseDiagnostic) && isscalar(D.Pointing.ChooseAntiSunNormal));
end

function X = validateHistoryVectors(X,nTimes,name)
% Explicit row convention avoids ambiguity when exactly three times are given.
    validateattributes(X,{'numeric'},{'2d','nonempty','real','finite'},mfilename,name);
    if nTimes==1 && isvector(X) && numel(X)==3, X=X(:)'; end
    assert(isequal(size(X),[nTimes 3]), ...
        '%s must be N-by-3: one Cartesian vector per time. Transpose 3-by-N arrays.',name);
end

function [Orbit,Q,r,v] = spacecraftFromCartesian(r,v,Earth)
% Use CDH Cartesian states directly, including radial velocity for eccentric
% spacecraft orbits. No element conversion, propagation, or SSO reset.
    radius = norm(r);
    assert(radius>Earth.Radius_km,'Spacecraft state is inside the Earth sphere.');
    angularMomentum = cross(r,v);
    momentum = norm(angularMomentum);
    assert(momentum>1e-12*radius*max(norm(v),1), ...
        'Cannot define RTN for a zero or nearly radial spacecraft velocity.');
    R = r/radius;
    N = angularMomentum/momentum;
    T = cross(N,R);
    Q = [R' T' N'];
    Orbit.Radius_km = radius;
    Orbit.Inclination_deg = acosd(max(-1,min(1,N(3))));
    % Normal velocity is zero by definition of this instantaneous orbit frame.
    Orbit.VelocityRTN_km_s = [dot(v,R),momentum/radius,0];
end

function stamp = fileStamp(path)
    assert(isfile(path),'Input file not found: %s',path);
    info = dir(path);
    assert(numel(info)==1,'Ambiguous input file: %s',path);
    stamp = {fullfile(info.folder,info.name),info.bytes,info.datenum};
end

function path = locateMASTERFile(Debris,Shell)
% Only filename selection is cached implicitly through file metadata later.
% Blank filename enables header-based discovery; normally use the pinned
% broad-range filename so older exported cases may remain in the folder.
    assert(isfolder(Debris.MASTERFolder),'MASTER folder not found: %s',Debris.MASTERFolder);
    if strlength(Debris.MASTERFileName)>0
        files = dir(fullfile(Debris.MASTERFolder,'**',char(Debris.MASTERFileName)));
        files = files(~[files.isdir]);
        assert(numel(files)==1,'Expected one %s beneath MASTERFolder; found %d.', ...
            Debris.MASTERFileName,numel(files));
        path = fullfile(files(1).folder,files(1).name);
    else
        files = dir(fullfile(Debris.MASTERFolder,'**','*'));
        matches = {};
        for k = 1:numel(files)
            if files(k).isdir || isempty(regexp(files(k).name,'\.__\d+$','once')), continue; end
            candidate = fullfile(files(k).folder,files(k).name);
            raw = fileread(candidate);
            x = regexp(raw,'x Data:\s*([^\r\n]+)','tokens','once');
            y = regexp(raw,'y Data:\s*([^\r\n]+)','tokens','once');
            if isempty(x) || isempty(y), continue; end
            if ~all(ismember({'S.D. Altitude [km]','S.D. Diameter [m]'}, ...
                    {strtrim(x{1}),strtrim(y{1})})), continue; end
            c = parseMASTERAltitudeDiameter(raw,candidate);
            if c.AltitudeEdges_km(1)<=Shell.MinAltitude_km && ...
                    c.AltitudeEdges_km(end)>=Shell.MaxAltitude_km && ...
                    c.DiameterEdges_cm(1)<=min(Debris.SizeBins_cm(:,1))+1e-7 && ...
                    c.DiameterEdges_cm(end)>=max(Debris.SizeBins_cm(:,2))-1e-7
                matches{end+1} = candidate; %#ok<AGROW>
            end
        end
        assert(numel(matches)==1,'Found %d candidate tables; set Debris.MASTERFileName explicitly.',numel(matches));
        path = matches{1};
    end
end

%% 6. SENSOR AND MASTER IMPORT HELPERS

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
    nBins = size(Debris.SizeBins_cm,1);
    assert(isfolder(Debris.MASTERFolder),'MASTER input folder does not exist: %s',Debris.MASTERFolder);
    path = fullfile(Debris.MASTERFolder,char(Debris.MASTERFileName));
    c = parseMASTERAltitudeDiameter(fileread(path),path);
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
    ih = discretize(altitude,Environment.AltitudeEdges_km);
    assert(all(isfinite(ih)),'Requested altitude is outside MASTER grid.');
    density = Environment.Density_per_km3(ih,bin);
end


%% 7. LOCAL ENCOUNTER GEOMETRY HELPERS

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
        for ih = 1:numel(edges)-1
            h = (edges(ih)+edges(ih+1))/2;
            debrisSpeed = sqrt(Earth.Mu_km3_s2/(Earth.Radius_km+h));
            for j = 1:numel(weights)
                phi = directions(j,1); node = directions(j,2);
                velocity = debrisSpeed*[0 cosd(phi) node*sind(phi)]-Orbit.VelocityRTN_km_s;
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
    end
end

function b = getBoresightRTN(Pointing,normalSign)
% RTN: R is radial out, T is along-track, N = R cross T.
    switch lower(Pointing.Mode)
        case "zenith", b = [1 0 0];
        case "alongtrack", b = [0 1 0];
        case "velocity", b = Pointing.CurrentVelocityRTN_km_s;
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

