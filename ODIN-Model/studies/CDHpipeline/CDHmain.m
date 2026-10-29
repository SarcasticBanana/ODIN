clear; close all; clc;

addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
odinPath;

cfg = config();
cfg.h = 600;                  % Your chosen altitude [km]
nOrbit = 5;
out = propagate(cfg, nOrbit);     % Ale's orbit calculation
N = length(out.t);

% lat = sign(out.r(3,:)) .* acosd(vecnorm(out.r(1:2,:)) ./ vecnorm(out.r));
% figure; plot(out.t,lat);
% figure; plot(lat,normpdf(lat,0,120));

[a_km, eccentricity, inclination_deg, RAAN_deg, argPerigee_deg, trueAnomaly_deg] = rv2oe(cfg,out);

%%
for i=1:N/10
    elements = [a_km(i), eccentricity(i), inclination_deg(i), RAAN_deg(i), argPerigee_deg(i), trueAnomaly_deg(i)];
    [rate_per_s,details] = debrisDetectionRate(elements,cfg,out.rsun(:,i));
    debrisRate(i) = details.TotalRate_per_s;
    disp(['Finished ', num2str(i)]);
end

debrisRate = repmat(debrisRate,1,10);
debrisRate = [debrisRate, debrisRate(end-3:end)];

%%
dataPerDebris = 10e6; % 10 Mb, guess
DLIntervals = [out.T*(0:nOrbit);out.T*(0:nOrbit)+out.T/10];
K_DL = 56e6; % 56 Mbps, probably
R_GNC = 100; % 100 b/s, guess
R_EPS = 100;
R_Struct = 10; % 10 b/s, guess
headerBits = 48;
dataBitsPerHeader = 65536*8;
packetizationInflation = (headerBits+dataBitsPerHeader)/dataBitsPerHeader;
pi_c = 5*packetizationInflation; % compression ratio, probably


opts = odeset('NonNegative',1,'MaxStep',10);
[tout,Dout] = ode45(@(t,D) dataTime(t,D,out.t,debrisRate,DLIntervals,dataPerDebris,K_DL,R_GNC,R_EPS,R_Struct,pi_c), ...
    [0, out.t(end)],[0; 0],opts);
%%
debrisRateOut = interp1(out.t,debrisRate,tout); % events/second

DLIdx = any(cell2mat(arrayfun(@(i) tout>=DLIntervals(1,i) & tout<=DLIntervals(2,i),1:length(DLIntervals),'UniformOutput',false))');

R_payload = debrisRateOut' * dataPerDebris .* double(~DLIdx);

% only downlink data that you can downlink
K_DLmin = min(Dout(:,1) + 1/pi_c*(R_GNC+R_EPS+R_Struct),K_DL)' .* double(DLIdx); 

dDdt = 1/pi_c*(R_payload+R_GNC+R_EPS+R_Struct) - K_DLmin;

figure; hold on;
plot(tout,Dout(:,1),'b','LineWidth',2);
ylabel("Data Stored (Bits)");
yyaxis('right');
plot(tout,dDdt,'r','LineWidth',2);
ylabel("Total Data Rate (Bits/s)");
xlabel("Time (s)");
box on; grid on;

function dDdt = dataTime(t,D,oldTime,debrisRate,DLIntervals,dataPerDebris,K_DL,R_GNC,R_EPS,R_Struct,pi_c)
    DLIdx = any(arrayfun(@(i) t>=DLIntervals(1,i) & t<=DLIntervals(2,i),1:length(DLIntervals)));

    R_payload = interp1(oldTime,debrisRate,t) * dataPerDebris * double(~DLIdx);

    % only downlink data that you can downlink
    K_DLmin = min(D(1) + 1/pi_c*(R_GNC+R_EPS+R_Struct),K_DL) * double(DLIdx);

    % disp(num2str(t));

    % if t>5797
    %     keyboard;
    % end
    dDdt = [1/pi_c*(R_payload+R_GNC+R_EPS+R_Struct) - K_DLmin;
            0];
end
