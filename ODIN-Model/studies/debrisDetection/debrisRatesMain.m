clear; clc; close all;

odinPath;

% read in ESA file and extract only data lines
raw = readlines("\studies\debrisDetection\MASTER_Data\altdecl_ODIN_2024_D010_050_t.txt");
fullRaw = raw(strlength(raw)~=0);
rawData = fullRaw(extract(fullRaw,1)~="#");

% extract altitudes, declinations, and densities
datalines = arrayfun(@(x) double(strsplit(x)),rawData,'UniformOutput',false);
dataMat = reshape([datalines{:}],[21 1440])';
alts = dataMat(:,2);
decs = deg2rad(dataMat(:,3)); % convert to radians!
dens = dataMat(:,end);

figure;
scatter3(alts,rad2deg(decs),dens,"filled","CData",dens);
ylim([-90 90]);
zscale('log');
grid on;
xlabel("Altitude Shell (km)");
ylabel("Declination (deg)");
zlabel("Spatial Density (1/km^3)");
title("Spatial Density at Declination and Altitude");

const.rE = 6378.137;

% calculate number of debris at each declination/altitude based on volume and density
rhos = const.rE+alts;
phis = pi/2-decs;
vols = 40*pi/3*sin(pi/72) * sin(phis).*(3*rhos.^2+25);
debrisCounts = vols.*dens;

figure;
scatter3(alts,rad2deg(decs),debrisCounts,"filled","CData",debrisCounts);
ylim([-90 90]);
grid on;
xlabel("Altitude Shell (km)");
ylabel("Declination (deg)");
zlabel("Number of Debris");
title("Debris Count at Declination and Altitude");
disp(['Counted ',num2str(round(sum(debrisCounts))),' total debris']);

%% reformat data into 2d matrix and plot over altitudes
close all;
altUniq = unique(alts);
decUniq = unique(decs);
nAlts = numel(altUniq);
nDecs = numel(decUniq);
debrisCounts2d = reshape(debrisCounts,[nDecs, nAlts]); % rows are declinations, cols are altitudes
debrisDens2d = reshape(dens,[nDecs, nAlts]); % rows are declinations, cols are altitudes

% figure;
% countPlt = plot(rad2deg(decUniq),debrisCounts2d(:,1),'b');
% ylim([0 max(debrisCounts)]);
% ylabel("Debris Counts");
% 
% yyaxis right;
% densPlt = plot(rad2deg(decUniq),debrisDens2d(:,1),'m');
% ylim([0 max(dens)]);
% ylabel("Spatial Density (1/km^3)");
% 
% xlim([-90 90]);
% xlabel("Declination (deg)");
% 
% ax = gca;
% ax.YAxis(1).Color = [0 0 1];
% ax.YAxis(2).Color = [1 0 1];
% 
% for i=2:nAlts
%     countPlt.YData = debrisCounts2d(:,i);
%     densPlt.YData = debrisDens2d(:,i);
%     title(sprintf("Debris Frequency at %i km",altUniq(i)));
%     pause(0.25);
%     drawnow limitrate
% end

%% calculate inclination distribution
nIncs = 500;
decs = decUniq(decUniq>0);
incs = linspace(0,pi,nIncs)';
[decsmesh,altsmesh,incsmesh] = meshgrid(decs,altUniq,incs);

Kmask = sin(incsmesh)>abs(sin(decsmesh));
K = zeros(nAlts,nDecs/2,nIncs);
K(Kmask) = 1./sqrt(max(1e-4,sin(incsmesh(Kmask)).^2-sin(decsmesh(Kmask)).^2));
dinc = median(diff(incs));

rhoModel = @(theta) densModel(theta,incs,dinc,K);

% Suppose rhoTarget is [numel(altUniq) x numel(posDecs)] on the same grid
rhoTarget = debrisDens2d(nDecs/2+1:end,:)';
validMask = ~isnan(rhoTarget) & isfinite(rhoTarget);
targetVec = rhoTarget(validMask);

mu0  = deg2rad([65, 74, 83, 98]);
A0   = ones(size(mu0)) * 1e-9;
sig0 = ones(size(mu0)) * deg2rad(3);
c0_init = 1e-12;

theta0 = [c0_init; reshape([A0; mu0; sig0], [], 1)];

lb = [0;     reshape([zeros(size(mu0)); zeros(size(mu0)); ones(size(mu0))*deg2rad(0.5)], [], 1)];
ub = [Inf;   reshape([Inf(size(mu0));   pi*ones(size(mu0));  ones(size(mu0))*deg2rad(30)], [], 1)];

resFun = @(th) localResiduals(th, incs, dinc, K, validMask, targetVec);

opts = optimoptions('lsqnonlin', 'Display', 'iter');
thetaHat = lsqnonlin(resFun, theta0, lb, ub, opts);

rhoHat = rhoModel(thetaHat);

function r = localResiduals(theta, incVec, dInc, K, validMask, targetVec)
    rhoPred = densModel(theta, incVec, dInc, K);
    r = rhoPred(validMask) - targetVec;
end

function rho = densModel(theta,incs,dinc,K)
    c0 = theta(1);
    params = reshape(theta(2:end), 3, []);   % rows: [A; mu; sig], cols: each lobe
    As   = params(1, :);
    mus  = params(2, :);
    sigmas = params(3, :);    
    winc = c0 + sum(As.*exp(-0.5*((incs-mus)./sigmas).^2),2);

    rho = sum(K .* reshape(winc, [1, 1, numel(incs)]), 3) * dinc;
end

%% plot inclination distribution for each altitude
close all;
figure;
for i=1:nAlts
    countsSing = debrisCounts2d(:,i);
    densSing = debrisDens2d(:,i);

    % decPDF = 1/(2*pi)*ones(1,nDecs); % use uniform for pdf
    % decPDF = countsSing/trapz(decUniq,countsSing); % use debris counts for pdf
    decPDF = densSing/trapz(decUniq,densSing); % use debris density for pdf

    % compute CDF for sampling
    decCDF = cumtrapz(decUniq,decPDF);
    N = 1e6;
    
    % sample declinations (pdf) and argument of latitudes (uniform)
    decSamples = interp1(decCDF, decUniq, rand(N,1), 'linear');
    randU = rand(N,1)*2*pi;
    
    % inclination of debris given declination at a current arugment of latitude
    % inc = asin(sin(dec)/sin(u))
    decOU = sin(decSamples)./sin(randU);
    decOU(decOU>1 | decOU<-1) = nan;
    decSamples(decOU>1 | decOU<-1) = nan;
    randU(decOU>1 | decOU<-1) = nan;
    incs = asin(decOU);
    
    % group inclinations into bins and calculate centers and counts
    % (to make sure aligns with decUniq)
    % histogram(incs,'BinWidth',5*pi/180);
    [incCounts(:,i), edges] = histcounts(incs,'BinWidth',5*pi/180);
    incCenters = edges(1:end-1)+2.5*pi/180;
    % pdf should just be the normalized counts distribution
    incPDF(:,i) = incCounts(:,i)/trapz(incCenters,incCounts(:,i));

    plot(rad2deg(incCenters),incPDF(:,i));
    xlim([-90 90]);
    ylim([0 0.5]);
    title(sprintf("Debris Frequency at %i km",altUniq(i)));
    xlabel("Inclination (deg)");
    ylabel("pdf");
    pause(0.5);
    drawnow limitrate;
end

%% observed declinations from spacecraft position (over time)
cfg = config();
cfg.h = 600;                  % Your chosen altitude [km]
nOrbit = 1;
out = propagate(cfg, nOrbit);     % Ale's orbit calculation
N = length(out.t);
% rsc = out.r(:,2);
% vsc = out.v(:,2);

% RECIs = reshape(out.r./vecnorm(out.r),[1,3,97]);
% NECIs = reshape(out.hvec./vecnorm(out.hvec),[1,3,97]);
% TECIs = cross(NECIs,RECIs);
% QECIs = [RECIs; TECIs; NECIs];
% 
% rsunsc = reshape(out.rsun-out.r,[3,1,97]);
% rsunRTN = pagemtimes(QECIs,rsunsc);

altsmesh = reshape(alts,[36,40]);
decsmesh = reshape(decs,[36,40]);

bRTN = [sind(66); 0; -cosd(66)];

alpha = deg2rad(1.49563646809545/2); % half width of sensor
xis = linspace(-1,1,10);
dxi = mean(diff(xis));
xis = xis(1:end-1)+dxi/2;
etas = linspace(-1,1,10);
deta = mean(diff(etas));
etas = etas(1:end-1)+deta/2;

rawSensorData = readtable("Sensor_Residence_Time_Trade 3.csv");
Rs = unique(rawSensorData.SlantRange_km);
dR = (circshift(Rs,-1)-circshift(Rs,1))/2;
dR(1) = dR(2);
dR(end) = dR(end-1);

xRTN = [0; 1; 0];
yRTN = cross(bRTN,xRTN);
QRTN = [xRTN yRTN bRTN];

decBins = decUniq-mean(diff(decUniq))/2;
decBins(end+1) = pi/2;

dV = zeros(length(xis),length(etas),length(Rs));
for it=1:length(out.t)
    disp(['looking at step ', num2str(it)])

    RhatECI = out.r(:,it)/norm(out.r(:,it));
    NECI = out.hvec(:,it)/norm(out.hvec(:,it));
    TECI = cross(NECI,RhatECI);
    QECI = [RhatECI TECI NECI];

    for i=1:length(xis)
        for j=1:length(etas)
            rcam(:,i,j) = [xis(i)*tan(alpha); etas(j)*tan(alpha); 1];
            usens(:,i,j) = rcam(:,i,j)/norm(rcam(:,i,j));
            uRTN(:,i,j) = QRTN*usens(:,i,j);
            for k=1:length(Rs)
                dV(i,j,k) = Rs(k)^2 * (tan(alpha)^2/(1+(xis(i)^2+etas(j)^2)*tan(alpha)^2)^2) * dR(k)*dxi*deta;
                RECI(:,i,j,k) = Rs(k)*QECI*uRTN(:,i,j);
                robsECI(:,i,j,k) = out.r(:,it) + RECI(:,i,j,k);
                hobsECI(i,j,k) = norm(robsECI(:,i,j,k)) - cfg.Re;
                if hobsECI(i,j,k)<605 || hobsECI(i,j,k)>995; hobsECI(i,j,k)=nan; end
                decobs(i,j,k) = asin(robsECI(3,i,j,k)/norm(robsECI(:,i,j,k)));
                densobs(i,j,k) = interp2(altsmesh,decsmesh,debrisDens2d,hobsECI(i,j,k),decobs(i,j,k));
                w(i,j,k) = densobs(i,j,k)*dV(i,j,k);
            end
        end
    end
    % disp(['min dec: ',num2str(min(decobs,[],"all"))]);
    % disp(['max dec: ',num2str(max(decobs,[],"all"))]);

    inShell(:,:,:,it) = (hobsECI>=605 & hobsECI<=995);
    
    decobs_valid = decobs(inShell(:,:,:,it));
    w_valid = w(inShell(:,:,:,it));
    binIdx = discretize(decobs_valid,decBins);
    validBin = ~isnan(binIdx);
    N_bin = accumarray(binIdx(validBin), w_valid(validBin), [numel(decBins)-1, 1], @sum, 0);

    p_dec_given_dsc(:,it) = N_bin ./ (sum(N_bin)*diff(decBins));

    plot(decUniq,p_dec_given_dsc(:,it));
    drawnow limitrate;
 
    % keyboard;
end

%% hella optimized
cfg = config();
cfg.h = 600;                  % Your chosen altitude [km]
nOrbit = 1;
out = propagate(cfg, nOrbit);     % Ale's orbit calculation
N = length(out.t);

altsmesh = reshape(alts,[36,40]);
decsmesh = reshape(decs,[36,40]);

bRTN = [sind(66); 0; -cosd(66)];

alpha = deg2rad(1.49563646809545/2); % half width of sensor
xis = linspace(-1,1,10);
dxi = mean(diff(xis));
xis = xis(1:end-1)+dxi/2;
etas = linspace(-1,1,10);
deta = mean(diff(etas));
etas = etas(1:end-1)+deta/2;

rawSensorData = readtable("Sensor_Residence_Time_Trade 3.csv");
Rs = unique(rawSensorData.SlantRange_km);
dR = (circshift(Rs,-1)-circshift(Rs,1))/2;
dR(1) = dR(2);
dR(end) = dR(end-1);

xRTN = [0; 1; 0];
yRTN = cross(bRTN,xRTN);
QRTN = [xRTN yRTN bRTN];

decBins = decUniq-mean(diff(decUniq))/2;
decBins(end+1) = pi/2;

% --- Precompute loop-invariant FOV geometry & differential volumes ---
nXi = numel(xis);
nEta = numel(etas);
nR = numel(Rs);

[xiGrid, etaGrid] = ndgrid(xis, etas);               % [nXi x nEta]
tanA = tan(alpha);
denom2D = sqrt(1 + (xiGrid.^2 + etaGrid.^2) * tanA^2);

% Unit ray vectors in camera frame and RTN frame: [3 x nXi x nEta]
usens = cat(3, xiGrid * tanA, etaGrid * tanA, ones(nXi, nEta)) ./ denom2D;
usens = permute(usens, [3, 1, 2]);
uRTN = pagemtimes(QRTN, usens);                      % [3 x nXi x nEta]

% Differential volume dV: [nXi x nEta x nR]
Rgrid3 = reshape(Rs, [1, 1, nR]);
Rgrid4 = reshape(Rs, [1, 1, 1, nR]);
angFactor2D = (tanA^2) ./ (denom2D.^4);
dV = angFactor2D .* (Rgrid3.^2) .* reshape(dR, [1, 1, nR]) * dxi * deta;

% Preallocate outputs across orbital steps
nDecBins = numel(decBins) - 1;
dDecBins = diff(decBins(:));
inShell = false(nXi, nEta, nR, N);
p_dec_given_dsc = zeros(nDecBins, N);

% Interpolant built once outside the time loop (transposed to ndgrid format)
densInterp = griddedInterpolant(altsmesh', decsmesh', debrisDens2d', 'linear', 'none');

% Flatten uRTN to [3 x nFOV] for fast matrix multiply per time step
nFOV = nXi * nEta;
uRTN_mat = reshape(uRTN, [3, nFOV]);

for it = 1:N
    disp(['looking at step ', num2str(it)])

    RhatECI = out.r(:,it) / norm(out.r(:,it));
    NECI = out.hvec(:,it) / norm(out.hvec(:,it));
    TECI = cross(NECI, RhatECI);
    QECI = [RhatECI TECI NECI];

    % Rotate all FOV unit rays into ECI and restore [3 x nXi x nEta] shape
    uECI = reshape(QECI * uRTN_mat, [3, nXi, nEta]);

    % Scale by slant range and add spacecraft position: [3 x nXi x nEta x nR]
    rRelECI = uECI .* Rgrid4;
    robsECI = out.r(:,it) + rRelECI;

    % Geocentric norm, altitude, and declination: [nXi x nEta x nR]
    robsNorm = reshape(vecnorm(robsECI, 2, 1), [nXi, nEta, nR]);
    hobsECI = robsNorm - cfg.Re;
    sinDec = max(-1, min(1, reshape(robsECI(3,:,:,:), [nXi, nEta, nR]) ./ robsNorm));
    decobs = asin(sinDec);

    % Mask to altitude shell before interpolating density
    mask = (hobsECI >= 605) & (hobsECI <= 995);
    inShell(:,:,:,it) = mask;

    w = zeros(nXi, nEta, nR);
    if any(mask, 'all')
        densValid = densInterp(hobsECI(mask), decobs(mask));
        validInterp = ~isnan(densValid);
        maskIdx = find(mask);
        w(maskIdx(validInterp)) = densValid(validInterp) .* dV(maskIdx(validInterp));
    end

    decobs_valid = decobs(mask);
    w_valid = w(mask);

    binIdx = discretize(decobs_valid, decBins);
    validBin = ~isnan(binIdx);
    N_bin = accumarray(binIdx(validBin), w_valid(validBin), [nDecBins, 1], @sum, 0);

    totalN = sum(N_bin);
    if totalN > 0
        p_dec_given_dsc(:,it) = N_bin ./ (totalN * dDecBins);
    end

    plot(decUniq, p_dec_given_dsc(:,it));
    drawnow limitrate;

    alphaobs = atan2(robsECI(2,:,:,:),robsECI(1,:,:,:));
    ehatr = robsECI ./ reshape(robsNorm,[1,nXi,nEta,nR]);
    ehatphi = [-sin(alphaobs); cos(alphaobs); zeros(size(alphaobs))];
    ehattheta = cross(ehatphi,ehatr);

    vdeb = sqrt(cfg.mu./robsNorm);

    % keyboard;
end

%% try to fit to model using nonlinear least squares or something
% this code is a failure. maybe we can try to make it work in the future.
% figure;
plot(decUniq,decPDF);

upperDec = decUniq(decUniq>0);
upperDebrisPDF = decPDF(decUniq>0);
decMesh = repmat(decUniq,1,nDecs);
incMesh = decMesh';
pdgiveni = 2*cos(decMesh)./(pi*sqrt(sin(incMesh).^2-sin(decMesh).^2));
%%
pdgiveni(decMesh>incMesh) = 0;
binWidth = mean(diff(decUniq));
pdgiveni(1:(length(pdgiveni)+1):end) = real(2/pi*(asin(sin(decMesh(:,1)+binWidth/2)./sin(incMesh(1,:)'))-asin(sin(decMesh(:,1)-binWidth/2)./sin(incMesh(1,:)'))));
%%
K_reg = pdgiveni; 
d_reg = upperDebrisPDF; 

% 4. Define the hard equality constraint: Area = sum(i_est * di) = 1
Aeq = ones(1, length(pdgiveni)) * binWidth;
beq = 1;

% 5. Define the lower bound (strictly non-negative probabilities)
lb = zeros(length(pdgiveni), 1);

% 6. Solve using lsqlin (suppressing the command window output)
options = optimoptions('lsqlin', 'Display', 'off');
i_est = lsqlin(K_reg, d_reg, [], [], Aeq, beq, lb, [], [], options);

figure; plot(upperDec,i_est);
