clear; close all; clc;

addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
odinPath;

cfg = config();
cfg.h = 600;                  % Your chosen altitude [km]
nOrbit = 3;
out = propagate(cfg, nOrbit);     % Ale's orbit calculation
N = length(out.t);

% lat = sign(out.r(3,:)) .* acosd(vecnorm(out.r(1:2,:)) ./ vecnorm(out.r));
% figure; plot(out.t,lat);
% figure; plot(lat,normpdf(lat,0,120));

% [rate_per_s,details] = debrisDetectionRate(out.t,out.r',out.v',vecnorm(out.r)'+cfg.Re,cfg,out.rsun');
% debrisRate = ones(1,N)*details.TotalRate_per_s;
load('debrisRates.mat');
debrisT = xData;
debrisRates = yData;
debrisRates = debrisRates(out.t<=out.T);
debrisRate = repmat(debrisRates,1,3);

dataPerDebris = 7.675e6; % 10 Mb, guess

% DLIntervals = [out.T*(0:nOrbit);out.T*(0:nOrbit)+out.T/10];
st = stations();
[dl,c] = downlinkTimes(out,cfg,st);
DLIntervals = dl';

K_DL = 30e6; % 56 Mbps, probably
R_GNC = 100; % 100 b/s, guess
R_EPS = 100;
R_Struct = 10; % 10 b/s, guess
headerBits = 48;
dataBitsPerHeader = 65536*8;
packetizationInflation = (headerBits+dataBitsPerHeader)/dataBitsPerHeader;
pi_c = 7.675/5.6*packetizationInflation; % compression ratio, probably

opts = odeset('NonNegative',1,'MaxStep',10);
[tout,Dout] = ode45(@(t,D) dataTime(t,D,out.t,debrisRate,DLIntervals,dataPerDebris,K_DL,R_GNC,R_EPS,R_Struct,pi_c), ...
    [0, out.t(end)],[0; 0],opts);
%
debrisRateOut = interp1(out.t,debrisRate,tout); % events/second

DLIdx = any(cell2mat(arrayfun(@(i) tout>=DLIntervals(1,i) & tout<=DLIntervals(2,i),1:length(DLIntervals),'UniformOutput',false))');

% R_payload = debrisRateOut' * dataPerDebris .* double(~DLIdx);
R_payload = dataPerDebris .* double(~DLIdx);

% only downlink data that you can downlink
K_DLmin = min(Dout(:,1) + 1/pi_c*(R_GNC+R_EPS+R_Struct),K_DL)' .* double(DLIdx); 

dDdt = 1/pi_c*(R_payload+R_GNC+R_EPS+R_Struct) - K_DLmin;

figure; hold on;
plot(tout,Dout(:,1),'b','LineWidth',2);
ylabel("Data Stored (Bits)");
yyaxis('right');
plot(tout,dDdt,'r','LineWidth',2);
yline(0,'r--');
ylabel("Total Data Rate (Bits/s)");
xlabel("Time (s)");
box on; grid on;

function dDdt = dataTime(t,D,oldTime,debrisRate,DLIntervals,dataPerDebris,K_DL,R_GNC,R_EPS,R_Struct,pi_c)
    DLIdx = any(arrayfun(@(i) t>=DLIntervals(1,i) & t<=DLIntervals(2,i),1:length(DLIntervals)));

    % R_payload = interp1(oldTime,debrisRate,t) * dataPerDebris * double(~DLIdx);
    R_payload = dataPerDebris * double(~DLIdx);

    % only downlink data that you can downlink
    K_DLmin = min(D(1) + 1/pi_c*(R_GNC+R_EPS+R_Struct),K_DL) * double(DLIdx);

    % disp(num2str(t));

    % if t>5797
    %     keyboard;
    % end
    dDdt = [1/pi_c*(R_payload+R_GNC+R_EPS+R_Struct) - K_DLmin;
            0];
end
