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
decUniqUpper = decUniq(decUniq>0);
incs = linspace(0,pi,nIncs)';
[decsmesh,altsmesh,incsmesh] = meshgrid(decUniqUpper,altUniq,incs);

Kmask = sin(incsmesh)>abs(sin(decsmesh));
K = zeros(nAlts,numel(decUniqUpper),nIncs);
K(Kmask) = 1./sqrt(max(1e-4,sin(incsmesh(Kmask)).^2-sin(decsmesh(Kmask)).^2));
dinc = median(diff(incs));

rhoModel = @(theta) densModel2D(theta,altUniq,incs,dinc,K);

% Suppose rhoTarget is [numel(altUniq) x numel(decUniqUpper)] on the same grid
rhoTarget = debrisDens2d(nDecs/2+1:end,:)';
validMask = ~isnan(rhoTarget) & isfinite(rhoTarget);
targetVec = rhoTarget(validMask);

% Initial guesses for 4 2D lobes: [A; mu_i; sig_i; mu_h; sig_h]
mu_i0  = deg2rad([65, 74, 83, 98]);
sig_i0 = ones(size(mu_i0)) * deg2rad(3);
mu_h0  = [800, 800, 830, 800];        % [km] anchored near the ~800 km peak
sig_h0 = ones(size(mu_i0)) * 50;      % [km]
A0     = ones(size(mu_i0)) * 1e-9;
c0_init = 1e-12;

theta0 = [c0_init; reshape([A0; mu_i0; sig_i0; mu_h0; sig_h0], [], 1)];

lb = [0;   reshape([zeros(size(mu_i0)); ...
                    zeros(size(mu_i0)); ...
                    ones(size(mu_i0))*deg2rad(0.5); ...
                    min(altUniq)*ones(size(mu_i0)); ...
                    10*ones(size(mu_i0))], [], 1)];

ub = [Inf; reshape([Inf(size(mu_i0)); ...
                    pi*ones(size(mu_i0)); ...
                    ones(size(mu_i0))*deg2rad(30); ...
                    max(altUniq)*ones(size(mu_i0)); ...
                    500*ones(size(mu_i0))], [], 1)];

resFun = @(th) localResiduals(th, altUniq, incs, dinc, K, validMask, targetVec);

opts = optimoptions('lsqnonlin', 'Display', 'iter');
thetaHat = lsqnonlin(resFun, theta0, lb, ub, opts);

rhoHat = rhoModel(thetaHat);

altsmesh2d = altsmesh(:,:,1);
decsmesh2d = decsmesh(:,:,1);
figure; hold on;
scatter3(altsmesh2d(:),decsmesh2d(:),rhoTarget(:),'b');
scatter3(altsmesh2d(:),decsmesh2d(:),rhoHat(:),'r');

function r = localResiduals(theta, altVec, incVec, dInc, K, validMask, targetVec)
    rhoPred = densModel2D(theta, altVec, incVec, dInc, K);
    r = 1e9*(rhoPred(validMask) - targetVec);
end

function winc2D = inclinationWeight(theta, altVec, incs)
    c0 = theta(1);
    params = reshape(theta(2:end), 5, []);   % rows: [A; mu_i; sig_i; mu_h; sig_h]
    As     = params(1, :);
    mu_i   = params(2, :);
    sig_i  = params(3, :);
    mu_h   = params(4, :);
    sig_h  = params(5, :);

    % Separable 2D Gaussians: G_h is [nAlts x nLobes], G_i is [nIncs x nLobes]
    G_h = exp(-0.5 * ((altVec(:) - mu_h) ./ sig_h).^2);
    G_i = exp(-0.5 * ((incs(:)   - mu_i) ./ sig_i).^2);

    % Combine into [nAlts x nIncs], then broadcast across declination (dim 2)
    winc2D = c0 + G_h * (As .* G_i)';
end

function rho = densModel2D(theta, altVec, incs, dinc, K)
    winc2D = inclinationWeight(theta, altVec, incs);
    rho = sum(K .* reshape(winc2D, [numel(altVec), 1, numel(incs)]), 3) * dinc;
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
out = propagate(cfg, nOrbit); % Ale's orbit calculation
N = length(out.t);

altsmesh = reshape(alts, [36, 40]);
decsmesh = reshape(decs, [36, 40]);

bRTN = [sind(66); 0; -cosd(66)];

alpha = deg2rad(1.49563646809545/2); % half width of sensor
xis = linspace(-1, 1, 10);
dxi = mean(diff(xis));
xis = xis(1:end-1) + dxi/2;
etas = linspace(-1, 1, 10);
deta = mean(diff(etas));
etas = etas(1:end-1) + deta/2;

rawSensorData = readtable("Sensor_Residence_Time_Trade 3.csv");
Rs = unique(rawSensorData.SlantRange_km);
dR = (circshift(Rs, -1) - circshift(Rs, 1)) / 2;
dR(1) = dR(2);
dR(end) = dR(end-1);

xRTN = [0; 1; 0];
yRTN = cross(bRTN, xRTN);
QRTN = [xRTN yRTN bRTN];

% --- Precompute loop-invariant FOV geometry & differential volumes ---
nXi  = numel(xis);
nEta = numel(etas);
nR   = numel(Rs);

[xiGrid, etaGrid] = ndgrid(xis, etas);               % [nXi x nEta]
tanA = tan(alpha);
denom2D = sqrt(1 + (xiGrid.^2 + etaGrid.^2) * tanA^2);

% Unit ray vectors in camera frame and RTN frame: [3 x nXi x nEta]
usens = cat(3, xiGrid * tanA, etaGrid * tanA, ones(nXi, nEta)) ./ denom2D;
usens = permute(usens, [3, 1, 2]);
uRTN  = pagemtimes(QRTN, usens);
uRTN_mat = reshape(uRTN, [3, nXi * nEta]);

% Differential volume dV: [nXi x nEta x nR]
Rgrid3 = reshape(Rs, [1, 1, nR]);
Rgrid4 = reshape(Rs, [1, 1, 1, nR]);
angFactor2D = (tanA^2) ./ (denom2D.^4);
dV = angFactor2D .* (Rgrid3.^2) .* reshape(dR, [1, 1, nR]) * dxi * deta;

% Interpolant built once outside the time loop (transposed to ndgrid format)
densInterp = griddedInterpolant(altsmesh', decsmesh', debrisDens2d', 'linear', 'none');

% Unpack fitted 2D inclination lobe parameters
c0     = thetaHat(1);
params = reshape(thetaHat(2:end), 5, []);
As     = params(1, :);
mu_i   = params(2, :);
sig_i  = params(3, :);
mu_h   = params(4, :);
sig_h  = params(5, :);

% Precompute inclination-side terms once outside the orbit loop
G_i     = exp(-0.5 * ((incs(:) - mu_i) ./ sig_i).^2); % [nIncs x nLobes]
incs4   = reshape(incs, [1, 1, 1, nIncs]);            % [1 x 1 x 1 x nIncs]
sinInc4 = sin(incs4);                                 % [1 x 1 x 1 x nIncs]
cosInc4 = cos(incs4);                                 % [1 x 1 x 1 x nIncs]

for it = 1:N
    disp(['looking at step ', num2str(it)])

    RhatECI = out.r(:,it) / norm(out.r(:,it));
    NECI    = out.hvec(:,it) / norm(out.hvec(:,it));
    TECI    = cross(NECI, RhatECI);
    QECI    = [RhatECI TECI NECI];

    % Rotate all FOV unit rays into ECI and scale by slant range: [3 x nXi x nEta x nR]
    uECI    = reshape(QECI * uRTN_mat, [3, nXi, nEta]);
    robsECI = out.r(:,it) + uECI .* Rgrid4;

    % Geocentric norm, altitude, and declination: [nXi x nEta x nR]
    robsNorm = reshape(vecnorm(robsECI, 2, 1), [nXi, nEta, nR]);
    hobsECI  = robsNorm - cfg.Re;
    sinDec   = max(-1, min(1, reshape(robsECI(3,:,:,:), [nXi, nEta, nR]) ./ robsNorm));
    decobs   = asin(sinDec);

    % Volumetric cell weight w: [nXi x nEta x nR] (zero outside shell)
    dens = densInterp(hobsECI, decobs);
    dens(isnan(dens) | hobsECI < 605 | hobsECI > 995) = 0;
    w = dens .* dV;

    % Local orthonormal basis at each cell: [3 x nXi x nEta x nR]
    alphaobs  = atan2(robsECI(2,:,:,:), robsECI(1,:,:,:));
    ehatr     = robsECI ./ reshape(robsNorm, [1, nXi, nEta, nR]);
    ehatphi   = [-sin(alphaobs); cos(alphaobs); zeros(size(alphaobs))];
    ehattheta = cross(ehatphi, ehatr, 1);

    % Circular debris speed and horizontal components: [nXi x nEta x nR x nIncs]
    vdeb   = sqrt(cfg.mu ./ robsNorm);
    ratio  = cosInc4 ./ cos(decobs);
    veast  = vdeb .* ratio;
    vnorth = vdeb .* sqrt(max(0, 1 - ratio.^2));

    % Debris and relative velocity in ECI (both branches): [3 x nXi x nEta x nR x nIncs]
    veast5  = reshape(veast,  [1, nXi, nEta, nR, nIncs]);
    vnorth5 = reshape(vnorth, [1, nXi, nEta, nR, nIncs]);

    vdebECI_asc  = veast5 .* ehatphi + vnorth5 .* ehattheta;
    vdebECI_desc = veast5 .* ehatphi - vnorth5 .* ehattheta;

    vsc = out.v(:,it);
    vrelECI_asc  = vdebECI_asc  - vsc;
    vrelECI_desc = vdebECI_desc - vsc;

    % Full N-D conditional inclination distribution p(i | h, dec): [nXi x nEta x nR x nIncs]
    G_h  = exp(-0.5 * ((hobsECI(:) - mu_h) ./ sig_h).^2);          % [(nXi*nEta*nR) x nLobes]
    winc = reshape(c0 + G_h * (As .* G_i)', [nXi, nEta, nR, nIncs]);

    diffSin2 = sinInc4.^2 - sin(decobs).^2;
    K_full   = zeros(size(diffSin2));
    Kmask    = sinInc4 > abs(sin(decobs));
    K_full(Kmask) = 1 ./ sqrt(max(1e-4, diffSin2(Kmask)));

    q   = K_full .* winc;
    p_i = q ./ sum(q, 4);
    p_i(isnan(p_i)) = 0;

    keyboard;
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
