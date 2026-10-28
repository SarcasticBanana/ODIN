clear; close all; clc;

addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
odinPath;

cfg = config();
cfg.h = 600;                  % Your chosen altitude [km]
nOrbit = 2;
out = propagate(cfg, nOrbit);     % Ale's orbit calculation
N = length(out.t);

lat = sign(out.r(3,:)) .* acosd(vecnorm(out.r(1:2,:)) ./ vecnorm(out.r));
% figure; plot(out.t,lat);

debrisRateCoeff = 1; % guess
dataPerDebris = 1e6; % 1 Mb, guess
DLIntervals = [out.T*(0:nOrbit);out.T*(0:nOrbit)+out.T/10];
K_DL = 1e4; % 1 Mbps, guess
R_GNC = 100; % 100 b/s, guess
R_EPS = 100;
R_Struct = 10; % 10 b/s, guess
pi_c = 5; % compression ratio, guess

opts = odeset('NonNegative',1,'MaxStep',10);
[tout,Dout] = ode45(@(t,D) dataTime(t,D,debrisRateCoeff,out.t,lat,DLIntervals,dataPerDebris,K_DL,R_GNC,R_EPS,R_Struct,pi_c), ...
    [0, out.t(end)],[0; 0],opts);

debrisRates = debrisRateCoeff*normpdf(interp1(out.t,lat,tout)',0,90); % events/second

DLIdx = any(cell2mat(arrayfun(@(i) tout>=DLIntervals(1,i) & tout<=DLIntervals(2,i),1:length(DLIntervals),'UniformOutput',false))');

R_payload = debrisRates * dataPerDebris .* double(~DLIdx);

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

function y = normpdf(x,mu,sigma)
    y = 1/(sigma*sqrt(2*pi)) * exp(-(x-mu).^2/(2*sigma^2));
end

function dDdt = dataTime(t,D,debrisRateCoeff,oldTime,lat,DLIntervals,dataPerDebris,K_DL,R_GNC,R_EPS,R_Struct,pi_c)
    debrisRates = debrisRateCoeff*normpdf(interp1(oldTime,lat,t),0,90); % events/second

    DLIdx = any(arrayfun(@(i) t>=DLIntervals(1,i) & t<=DLIntervals(2,i),1:length(DLIntervals)));

    R_payload = debrisRates * dataPerDebris * double(~DLIdx);

    % only downlink data that you can downlink
    K_DLmin = min(D(1) + 1/pi_c*(R_GNC+R_EPS+R_Struct),K_DL) * double(DLIdx);

    % disp(num2str(t));

    % if t>5797
    %     keyboard;
    % end
    dDdt = [1/pi_c*(R_payload+R_GNC+R_EPS+R_Struct) - K_DLmin;
            0];
end
