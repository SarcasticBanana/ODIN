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

    plot(rad2deg(incCentzers),incPDF(:,i));
    xlim([-90 90]);
    ylim([0 0.5]);
    title(sprintf("Debris Frequency at %i km",altUniq(i)));
    xlabel("Inclination (deg)");
    ylabel("pdf");
    pause(0.5);
    drawnow limitrate;
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
