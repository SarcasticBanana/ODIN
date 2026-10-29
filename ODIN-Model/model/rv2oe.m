function [a_km, eccentricity, inclination_deg, RAAN_deg, argPerigee_deg, trueAnomaly_deg] = rv2elems(cfg,out)
    % Using Algorithm 4.2 in Orbital Mechanics for Engineering Students
    % r: inertial position vector (XYZ, rows) over time (columns)
    % v: inertial velocity vector (XYZ, rows) over time (columns)
    dist = vecnorm(out.r);
    speed = vecnorm(out.v);
    radialVel = dot(out.r,out.v)./dist;
    htot = vecnorm(out.hvec);
    inclination_deg = acosd(out.hvec(3,:)./htot);

    N = cross([zeros(1,length(out.hvec));zeros(1,length(out.hvec));ones(1,length(out.hvec))],out.hvec);
    Ntot = vecnorm(N);
    try RAAN_deg(N(2,:)>=0) = acosd(N(1,N(2,:)>=0)./Ntot);
    catch, disp("Node vec Y is never positive"); end
    try RAAN_deg(N(2,:)<0) = 360 - acosd(N(1,N(2,:)>=0)./Ntot);
    catch, disp("Node vec Y is never negative"); end

    evec = 1/cfg.mu * ((speed.^2 - cfg.mu./dist).*out.r - dist.*radialVel.*out.r);
    % eccentricity = real(sqrt(1 + htot.^2/cfg.mu^2 .* (speed.^2 - 2*cfg.mu./dist)));
    eccentricity = vecnorm(evec);

    posEZ = evec(3,:)>=0;
    negEZ = evec(3,:)<0;
    try argPerigee_deg(posEZ) = acosd(dot(N(:,posEZ), evec(:,posEZ))./(Ntot(posEZ).*eccentricity(posEZ)));
    catch, disp("Eccentricity vec Z is never positive"); end
    try argPerigee_deg(negEZ) = 360 - acosd(dot(N(:,negEZ), evec(:,negEZ))./(Ntot(negEZ).*eccentricity(negEZ)));
    catch, disp("Eccentricity vec Z is never negative"); end

    posVR = radialVel>=0 & ~isapprox(eccentricity,0,'veryloose');
    negVR = radialVel<0 & ~isapprox(eccentricity,0,'veryloose');
    posRZ = out.r(3,:)>=0 & isapprox(eccentricity,0,'veryloose');
    negRZ = out.r(3,:)<=0 & isapprox(eccentricity,0,'veryloose');

    try trueAnomaly_deg(posVR) = acosd(1./eccentricity(posVR) .* (htot(posVR).^2./(cfg.mu*dist(posVR)) - 1));
    catch, disp("Radial velocity is never positive or eccentricity is never nonzero"); end
    try trueAnomaly_deg(negVR) = 360 - acosd(1./eccentricity(negVR) .* (htot(negVR).^2./(cfg.mu*dist(negVR)) - 1));
    catch, disp("Radial velocity is never negative or eccentricity is never nonzero"); end
    try trueAnomaly_deg(posRZ) = acosd(dot(N(:,posRZ),out.r(:,posRZ)) ./ (Ntot(posRZ).*dist(posRZ)));
    catch, disp("R vec Z is never negative or eccentricity is never zero"); end
    try trueAnomaly_deg(negRZ) = 360 - acosd(dot(N(:,negRZ),out.r(:,negRZ)) ./ (Ntot(negRZ).*dist(negRZ)));
    catch, disp("R vec Z is never negative or eccentricity is never zero"); end

    rperigee = htot.^2./(cfg.mu*(1+eccentricity));
    rapogee = htot.^2./(cfg.mu*(1-eccentricity));

    a_km = 0.5*(rperigee+rapogee);
end

