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
nIncs = 1000;
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
                    ones(size(mu_i0))*deg2rad(0.1); ...
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
view(3);

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
    incSamples = asin(decOU);
    
    % group inclinations into bins and calculate centers and counts
    % (to make sure aligns with decUniq)
    % histogram(incs,'BinWidth',5*pi/180);
    [incCounts(:,i), edges] = histcounts(incSamples,'BinWidth',5*pi/180);
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

%% observed relative velocity and slant ranges from spacecraft position (over time) (hella optimized)
cfg = config();
cfg.h = 600;
nOrbit = 1;
out = propagate(cfg, nOrbit);
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

sizeFile = "\studies\debrisDetection\MASTER_Data\altsize_ODIN_2024_D010_070_t.txt";
qualify = buildDebrisQualification(sizeFile, rawSensorData);

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

% Define velocity [km/s] bin edges
vBins  = 0:0.25:16;
dv     = diff(vBins(:));             % [nVBins x 1]
nVBins = numel(vBins) - 1;

% Preallocate outputs
p_vrel           = zeros(nVBins, nR, N);
NqualHistory     = zeros(nVBins, nR, N); % qualified objects per speed-range bin
NinView          = zeros(N, 1);          % all 1-5 cm objects before qualification
NqualInView      = zeros(N, 1);          % qualified objects in view (all speeds)
NqualOutsideBins = zeros(N, 1);          % qualified weight outside vBins limits

% Precompute slant-range bin index matching [nXi x nEta x nR x nIncs]
rIdx4 = repmat(reshape(1:nR, [1, 1, nR, 1]), [nXi, nEta, 1, nIncs]);

stopToc = 0;
totalToc = 0;

% Pre-allocate for diagnostic check
% Population-weighted means for all modeled 1–5 cm debris in view.
% NaN means the sample has no modeled population from which to form a mean.
meanRange_km   = nan(N,1);
meanVperp_km_s = nan(N,1);

for it = 1:N
    fprintf("looking at step %02i/%i (%i%%, %.1fs left)\n",it,N,round(100*it/N),totalToc/(it/N) - totalToc)
    tic;

    rSC = out.r(:,it);
    vsc = out.v(:,it);

    RhatECI = rSC / norm(rSC);
    NECI    = out.hvec(:,it) / norm(out.hvec(:,it));
    TECI    = cross(NECI, RhatECI);
    QECI    = [RhatECI TECI NECI];

    % Rotate all FOV unit rays into ECI and scale by slant range: [3 x nXi x nEta x nR]
    uECI    = reshape(QECI * uRTN_mat, [3, nXi, nEta]);
    robsECI = rSC + uECI .* Rgrid4;

    % Geocentric norm, altitude, and declination: [nXi x nEta x nR]
    robsNorm = reshape(vecnorm(robsECI, 2, 1), [nXi, nEta, nR]);
    hobsECI  = robsNorm - cfg.Re;
    sinDec   = max(-1, min(1, reshape(robsECI(3,:,:,:), [nXi, nEta, nR]) ./ robsNorm));
    decobs   = asin(sinDec);

    %Diagnostic check for debris-cell declinations actually reached
    inAltitude = hobsECI >= 605 & hobsECI <= 995;
    polarEdgeVolume_km3(it) = sum(dV(inAltitude & abs(decobs) > max(decUniq)));

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

    % Branch weights and scalar relative speeds: [nXi x nEta x nR x nIncs]
    w_branch     = 0.5 * w .* p_i;
    vrelMag_asc  = reshape(vecnorm(vrelECI_asc,  2, 1), [nXi, nEta, nR, nIncs]);
    vrelMag_desc = reshape(vecnorm(vrelECI_desc, 2, 1), [nXi, nEta, nR, nIncs]);

    % Angular velocity of camera fixed in RTN [rad/s] and frame-rotation velocity [km/s]
    omegaECI    = cross(rSC, vsc) / dot(rSC, rSC);
    omegaOnRays = repmat(omegaECI, [1, nXi, nEta]);
    vFrame      = cross(omegaOnRays, uECI, 1) .* Rgrid4;

    % Vectorized apparent sky-plane (transverse) speed [km/s]: [nXi x nEta x nR x nIncs]
    vAppAsc    = vrelECI_asc  - vFrame;
    vParAsc    = sum(uECI .* vAppAsc, 1);
    vPerp_asc  = reshape(sqrt(max(0, sum(vAppAsc.^2, 1)  - vParAsc.^2)),  [nXi, nEta, nR, nIncs]);

    vAppDesc   = vrelECI_desc - vFrame;
    vParDesc   = sum(uECI .* vAppDesc, 1);
    vPerp_desc = reshape(sqrt(max(0, sum(vAppDesc.^2, 1) - vParDesc.^2)), [nXi, nEta, nR, nIncs]);

    % Keep only nonzero-weight entries (already filtered to the 605-995 km shell)
    nz = w_branch > 0;
    if any(nz, 'all')
        vAll = [vrelMag_asc(nz); vrelMag_desc(nz)];
        rAll = [rIdx4(nz);       rIdx4(nz)];
        wAll = [w_branch(nz);    w_branch(nz)];

        % Map each retained sample back to its spatial cell for altitude
        branchIdx = find(nz);
        cellIdx   = mod(branchIdx - 1, numel(hobsECI)) + 1;
        altAll    = [hobsECI(cellIdx); hobsECI(cellIdx)];

        speedForSensor = [vPerp_asc(nz); vPerp_desc(nz)];
        pSizeQual      = qualify(altAll, rAll, speedForSensor);
        wQual          = wAll .* pSizeQual;

        % Expected objects in view before and after sensor qualification
        NinView(it)     = sum(wAll);
        NqualInView(it) = sum(wQual);

        % rAll contains range indices. Convert to physical slant range before averaging.
        meanRange_km(it)   = dot(wAll, Rs(rAll))      / NinView(it);
        meanVperp_km_s(it) = dot(wAll, speedForSensor) / NinView(it);

        vIdx = discretize(vAll, vBins);
        ok   = ~isnan(vIdx);

        % Accumulate across xi, eta, altitude, inclination, and asc/desc branches
        binSubs  = [vIdx(ok), rAll(ok)];
        N_vr     = accumarray(binSubs, wAll(ok),  [nVBins, nR], @sum, 0);
        Nqual_vr = accumarray(binSubs, wQual(ok), [nVBins, nR], @sum, 0);

        NqualHistory(:,:,it) = Nqual_vr;
        NqualOutsideBins(it) = sum(wQual(~ok));
        fracQual_vr = Nqual_vr ./ N_vr;
        fracQual_vr(N_vr == 0) = NaN;

        % Normalize along dim 1 (v_rel) to get p(v_rel | R, t)
        pdf_it = N_vr ./ (sum(N_vr, 1) .* dv);
        pdf_it(isnan(pdf_it)) = 0;
        p_vrel(:,:,it) = pdf_it;
    end
    stopToc = toc;
    totalToc = totalToc+stopToc;
end

%% Whole-orbit histogram and relative speed versus slant range animation
% Rerun this section for playback with the completed histories in memory.
framePause_s = 0.08;
nFrames = size(p_vrel,3);
t = out.t(:);
if isdatetime(t) || isduration(t)
    t = seconds(t-t(1));
else
    t = double(t-t(1));
end
assert(numel(t)==nFrames && all(isfinite(t)) && all(diff(t)>0), ...
    'Expected one increasing timestamp per completed frame.');

% Sum qualified occupancy over range; average over elapsed orbit time.
bySpeed = reshape(sum(NqualHistory,2),numel(vBins)-1,nFrames);
if nFrames == 1
    meanBySpeed = bySpeed;
else
    meanBySpeed = trapz(t,bySpeed,2)/(t(end)-t(1));
end

% Replace only earlier copies of these two outputs when replaying.
delete(findall(groot,'Type','figure','Tag','debrisOrbitHistogram'));
delete(findall(groot,'Type','figure','Tag','debrisOrbitAnimation'));
fHist = figure('Name','Whole-orbit qualified occupancy','Color','w','Tag','debrisOrbitHistogram');
ax = axes(fHist);
hold on;
% histogram(ax,'BinEdges',vBins500,'BinCounts',meanBySpeed500.','FaceAlpha',0.25);
% histogram(ax,'BinEdges',vBins1k,'BinCounts',meanBySpeed1k.','FaceAlpha',0.25);
% legend("nIncs = 500","nIncs = 1000");
histogram(ax,'BinEdges',vBins,'BinCounts',meanBySpeed.','FaceColor',[0.10 0.42 0.65]);
xlabel(ax,'Total relative speed (km/s)');
ylabel(ax,'Mean qualified objects in view per bin');
title(ax,'Qualified debris over the whole orbit','Color','k');
grid(ax,'on'); ax.Color='w'; ax.XColor='k'; ax.YColor='k';

% fAnim = figure('Name','Relative speed versus slant range','Color','w','Tag','debrisOrbitAnimation');
% ax = axes(fAnim);
% rangeEdges = [Rs(1)-(Rs(2)-Rs(1))/2; (Rs(1:end-1)+Rs(2:end))/2; Rs(end)+(Rs(end)-Rs(end-1))/2];
% % pcolor requires one extra vertex on each axis to show every bin.
% h = pcolor(ax,vBins,rangeEdges,p_vrel([1:end end],[1:end end],1).');
% h.EdgeColor = 'none';
% ax.YDir='normal'; ax.Color='w'; ax.XColor='k'; ax.YColor='k';
% clim(ax,[0 max(eps,max(p_vrel,[],'all'))]);
% cb = colorbar(ax); cb.Color='k'; cb.Label.String='Conditional speed density (s/km)';
% xlabel(ax,'Total relative speed (km/s)'); ylabel(ax,'Slant range (km)');
% heading = title(ax,'','Color','k');
% view(3);
% for frame = 1:nFrames
%     if ~isgraphics(h), break; end
%     h.ZData = p_vrel([1:end end],[1:end end],frame).';
%     h.CData = p_vrel([1:end end],[1:end end],frame).';
%     zlim([0 2]);
%     heading.String = sprintf('Relative speed at each range | Orbit sample %d of %d',frame,nFrames);
%     drawnow;
%     pause(framePause_s);
% end


%% Qualified objects in the FOV versus spacecraft argument of latitude
% Optional: replace [] with a logical vector aligned with out.t.
% Example: contactMask = groundStationInView(:);
contactMask = [];

% Recover the total from the speed-range bins at each orbit sample.
% NqualOutsideBins includes qualified objects outside the plotted speed range.
time = out.t(:);
nSamples = numel(time);
qualifiedFromBins = reshape(sum(NqualHistory, [1 2]), [], 1) ...
    + NqualOutsideBins(:);
qualifiedInFOV = NqualInView(:);

% The two independently accumulated totals should agree to roundoff.
maxDifference = max(abs(qualifiedFromBins - qualifiedInFOV));
fprintf('Maximum difference between direct and binned totals: %.3g objects\n', ...
    maxDifference);

% Argument of latitude u is measured in the direction of motion from the
% ascending node. It remains well defined for this inclined circular orbit.
% Each column of out.r and out.hvec corresponds to one time sample.
rHat = out.r ./ vecnorm(out.r, 2, 1);
hHat = out.hvec ./ vecnorm(out.hvec, 2, 1);

% Ascending-node direction = ECI north pole crossed with orbit normal.
nodeECI = [-hHat(2,:); hHat(1,:); zeros(1,nSamples)];
nodeHat = nodeECI ./ vecnorm(nodeECI, 2, 1);
quarterOrbitHat = cross(hHat, nodeHat, 1);

uDeg = mod(atan2d(sum(rHat .* quarterOrbitHat, 1), ...
    sum(rHat .* nodeHat, 1)), 360).';

% Keep the table in time order for comparison with communications windows.
orbitOpportunity = table((1:nSamples).', time, uDeg, qualifiedInFOV, ...
    'VariableNames', {'Sample','Time','ArgumentOfLatitude_deg','QualifiedInFOV'});

if ~isempty(contactMask)
    assert(islogical(contactMask) && numel(contactMask) == nSamples, ...
        'contactMask must contain one logical value per orbit sample.');
    orbitOpportunity.InContact = contactMask(:);
end

% Sort only the plotted coordinates so the line does not jump from 360 to 0.
% The current driver propagates one orbit (nOrbit = 1).
[uPlot, order] = sort(uDeg);
[peakQualified, peakIdx] = max(qualifiedInFOV);

figure('Name','Qualified objects versus orbit position','Color','w');
ax = axes;
plot(ax, uPlot, qualifiedInFOV(order), '-', 'LineWidth', 1.6, ...
    'Color', [0.10 0.36 0.65], 'DisplayName', 'Qualified objects in FOV');
hold(ax,'on');

if ~isempty(contactMask) && any(contactMask)
    plot(ax, uDeg(contactMask), qualifiedInFOV(contactMask), 'o', ...
        'Color', [0.85 0.35 0.08], 'MarkerSize', 5, ...
        'DisplayName', 'During communications contact');
end

plot(ax, uDeg(peakIdx), peakQualified, 'kp', ...
    'MarkerFaceColor', 'y', 'MarkerSize', 11, 'DisplayName', 'Sampled peak');
xlim(ax,[0 360]);
xticks(ax,0:45:360);
grid(ax,'on');
ax.Color = 'w';
ax.XColor = 'k';
ax.YColor = 'k';
xlabel(ax,'Spacecraft argument of latitude, u (deg)');
ylabel(ax,'Expected qualified objects in FOV');
title(ax,'Qualified debris in view versus orbit position','Color','k');
legend(ax,'Location','best');

fprintf('\nPeak sampled occupancy: %.6g qualified objects at u = %.2f deg\n', ...
    peakQualified, uDeg(peakIdx));
%disp(orbitOpportunity(peakIdx,:));

%% Diagnose the high-latitude dips
% Fraction of modeled in-view debris that passes Hunter's sensor threshold.
qualFraction = NqualInView(:) ./ NinView(:);
qualFraction(NinView(:) == 0) = NaN;

figure('Color','w','Name','Polar occupancy diagnostic');
tiledlayout(2,1,'TileSpacing','compact');

ax1 = nexttile;
plot(uDeg(order), NinView(order), '--', ...
    uDeg(order), NqualInView(order), '-', 'LineWidth',1.5);
ylabel('Expected objects in FOV');
legend('All modeled 1–5 cm debris','Sensor-qualified debris','Location','best');
title('Does the dip appear before sensor qualification?');

ax2 = nexttile;
plot(uDeg(order), 100*qualFraction(order), 'LineWidth',1.5);
ylabel('Qualified fraction (%)');
xlabel('Spacecraft argument of latitude (deg)');

for ax = [ax1 ax2]
    xline(ax,90,':','HandleVisibility','off');
    xline(ax,270,':','HandleVisibility','off');
    xlim(ax,[0 360]); xticks(ax,0:45:360); grid(ax,'on');
    ax.Color = 'w'; ax.XColor = 'k'; ax.YColor = 'k';
end

%% Compare population-weighted range and apparent speed by orbit position
% Requires uDeg from the argument-of-latitude section.
[uPlot, order] = sort(uDeg(:));
qualifiedPct = 100 * NqualInView(:) ./ NinView(:);
qualifiedPct(NinView(:) == 0) = NaN;

% Report actual sampled angles because the orbit grid may miss the targets.
targetDeg = [90; 100; 260; 270];
angleError = abs(mod(uDeg(:) - targetDeg.' + 180, 360) - 180);
[~, idx] = min(angleError, [], 1);
idx = idx(:);

polarComparison = table(targetDeg, uDeg(idx), meanRange_km(idx), ...
    meanVperp_km_s(idx), qualifiedPct(idx), ...
    'VariableNames', {'TargetU_deg','SampledU_deg','MeanRange_km', ...
    'MeanVperp_km_s','Qualified_pct'});
%disp(polarComparison);

figure('Color','w','Name','Range, transverse speed, and qualification');
layout = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
series = [meanRange_km, meanVperp_km_s, qualifiedPct];
labels = {'Mean slant range (km)', 'Mean transverse speed (km/s)', ...
    'Qualified fraction (%)'};

for j = 1:3
    ax = nexttile(layout);
    plot(ax, uPlot, series(order,j), 'LineWidth',1.5);
    hold(ax,'on');
    plot(ax, uDeg(idx), series(idx,j), 'ko', 'MarkerFaceColor','w');
    xline(ax,90,':');
    xline(ax,270,':');
    xlim(ax,[0 360]);
    xticks(ax,0:45:360);
    ylabel(ax,labels{j});
    grid(ax,'on');
    ax.Color = 'w'; ax.XColor = 'k'; ax.YColor = 'k';
end
xlabel(ax,'Spacecraft argument of latitude (deg)');
title(layout,'Debris in the FOV before sensor qualification','Color','k');


%% Downlink windows and missed qualified observing exposure

% Optional: supply a qualified passage rate [detections/s], aligned with out.t.
% Leave empty when only NinView and NqualInView are available.
qualifiedRate_per_s = [];

% Default station set from stations.m. Overlapping contacts are already merged.
[dl, stationContacts] = downlinkTimes(out, cfg);
dl = dl(dl(:,2) > dl(:,1), :);       % Keep positive-duration intervals
tSec = out.t(:);                     % Same seconds and time origin as dl
nQual = NqualInView(:);
assert(numel(tSec) == numel(nQual) && numel(tSec) >= 2 && all(diff(tSec) > 0), ...
    'Expected one qualified occupancy value per increasing orbit timestamp.');

% Unwrap in TIME order so interpolation across 359 -> 0 deg follows the orbit.
% The axis tick labels below are converted back to conventional 0–360 deg.
uContinuous = rad2deg(unwrap(deg2rad(uDeg(:))));

hasRate = ~isempty(qualifiedRate_per_s);
if hasRate
    qualifiedRate_per_s = qualifiedRate_per_s(:);
    assert(numel(qualifiedRate_per_s) == numel(tSec) && ...
        all(isfinite(qualifiedRate_per_s) & qualifiedRate_per_s >= 0), ...
        'Supply one finite, nonnegative passage rate per orbit timestamp.');
end

% Integrate to the actual contact boundaries, including short contacts
% that begin and end between two propagation samples.
nWindows = size(dl,1);
duration_s = dl(:,2) - dl(:,1);
lostExposure = zeros(nWindows,1);     % qualified object-seconds
lostDetections = zeros(nWindows,1);   % used only when a passage rate is supplied

for k = 1:nWindows
    aos = dl(k,1);
    los = dl(k,2);
    inside = tSec > aos & tSec < los;
    tw = [aos; tSec(inside); los];

    lostExposure(k) = trapz(tw, interp1(tSec, nQual, tw, 'linear'));
    if hasRate
        lostDetections(k) = trapz(tw, ...
            interp1(tSec, qualifiedRate_per_s, tw, 'linear'));
    end
end

% Report where each contact occurs and what it costs.
uStart = interp1(tSec, uContinuous, dl(:,1), 'linear');
uEnd   = interp1(tSec, uContinuous, dl(:,2), 'linear');
downlinkReport = table(dl(:,1), dl(:,2), duration_s, ...
    mod(uStart,360), mod(uEnd,360), lostExposure, ...
    'VariableNames', {'Start_s','End_s','Duration_s','StartU_deg', ...
    'EndU_deg','MissedQualifiedObjectSeconds'});
if hasRate
    downlinkReport.ExpectedMissedDetections = lostDetections;
end
disp(downlinkReport);

orbitExposure = trapz(tSec, nQual);
fprintf('\nDownlink time: %.2f min (%.2f%% of modeled time)\n', ...
    sum(duration_s)/60, 100*sum(duration_s)/(tSec(end)-tSec(1)));
fprintf('Qualified exposure missed: %.5g of %.5g object-s\n', ...
    sum(lostExposure), orbitExposure);
if orbitExposure > 0
    fprintf('Qualified exposure lost: %.2f%%\n', ...
        100*sum(lostExposure)/orbitExposure);
end
if hasRate
    fprintf('Expected missed detections: %.4g\n', sum(lostDetections));
end

% Plot the potential observing opportunity, including what is lost during contact.
figure('Color','w','Name','Qualified debris and downlink windows');
ax = axes;
hold(ax,'on');
yTop = 1.08*max(nQual);
if yTop <= 0, yTop = 1; end

for k = 1:nWindows
    band = patch(ax, [uStart(k) uEnd(k) uEnd(k) uStart(k)], ...
        [0 0 yTop yTop], [0.95 0.66 0.28], ...
        'FaceAlpha',0.30,'EdgeColor','none', ...
        'DisplayName','Downlink / payload off');
    if k > 1, band.HandleVisibility = 'off'; end
end

plot(ax,uContinuous,nQual,'LineWidth',1.7,'Color',[0.10 0.35 0.65], ...
    'DisplayName','Qualified objects in FOV');
xlim(ax,[uContinuous(1) uContinuous(end)]);
ylim(ax,[0 yTop]);
ticks = (45*ceil(uContinuous(1)/45)):45:(45*floor(uContinuous(end)/45));
xticks(ax,ticks);
xticklabels(ax,string(mod(ticks,360)));
xlabel(ax,'Spacecraft argument of latitude (deg)');
ylabel(ax,'Expected qualified objects in FOV');
title(ax,'Qualified debris opportunities during downlink','Color','k');
ax.Color = 'w'; ax.XColor = 'k'; ax.YColor = 'k';
grid(ax,'on'); legend(ax,'Location','best');