clear; close all; clc;

addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
odinPath;

cfg = config();
cfg.h = 600;                  % Your chosen altitude [km]
nOrbit = 10;
out = propagate(cfg, nOrbit);     % Ale's orbit calculation
N = length(out.t);

lat = sign(out.r(3,:)) .* acosd(vecnorm(out.r(1:2,:)) ./ vecnorm(out.r));

debrisRateCoeff = 1; % guess
dataPerDebris = 1e6; % 1 Mb, guess
%DLIntervals = [out.T*(0:nOrbit);out.T*(0:nOrbit)+out.T/10];
% K_DL = 1e4; % 1 Mbps, guess
% R_GNC = 100; % 100 b/s, guess
% R_EPS = 100;
% R_Struct = 10; % 10 b/s, guess
% pi_c = 5; % compression ratio, guess
% 
% opts = odeset('NonNegative',1,'MaxStep',10);
% [tout,Dout] = ode45(@(t,D) dataTime(t,D,debrisRateCoeff,out.t,lat,DLIntervals,dataPerDebris,K_DL,R_GNC,R_EPS,R_Struct,pi_c), ...
%     [0, out.t(end)],[0; 0],opts);

noiseSamples = load('model/Payload Models/Event Data Simulation Model/eventSamplesNoise.mat');
dataRateTimes = noiseSamples.tspost;
dataRates = noiseSamples.dataRate16;

% DRTimesInterp = linspace(0,max(dataRateTimes),100);
% figure; hold on;
% scatter(dataRateTimes,dataRates/1e6);
% plot(DRTimesInterp,interp1(dataRateTimes,dataRates/1e6,DRTimesInterp,'linear','extrap'),'-');
% xlabel('Time (s)'); ylabel('Data Rate (Mbps');
% grid on; box on;

st = stations();
[dl,c] = downlinkTimes(out,cfg,st);
DLIntervals = dl';

K_DL = 56e6; % 56 Mbps, probably
R_GNC = 230.4e3; % XACT baud: 230.4kHz
R_EPS = 100; % 100 b/s, guess
R_Struct = 100;
headerBits = 48;
dataBitsPerHeader = 65536*8;
packetizationInflation = (headerBits+dataBitsPerHeader)/dataBitsPerHeader;
pi_cAvg = 2.2357;
pi_c = @(x) (randn(1,x)*0.0282+pi_cAvg)/packetizationInflation; % compression ratio,
pi_cVal = pi_cAvg/packetizationInflation; %maube
storage_cap = 0.8*64 * 8e9; % [b] Upper bound on storage SD card
recollect_threshold = 0; % [b] Describes threshold required during calibration to start re-collecting data after SD card fills up.

dt = 0.5;
%global eventHappening tEvent downlinking
%eventHappening = 0;
%tEvent = 0;
%downlinking = 0;

phase_type = "Calibrate";
% if strcmpi(phase_type, "Calibrate")
%     opts = odeset('NonNegative',1,'MaxStep',10,'AbsTol',1); % just collect data, no need to model probabilist events
%     t_vec = [0, out.t(end)]; % Let ode45 choose step
% else
%     opts = odeset('NonNegative',1,'MaxStep',10,'AbsTol',1);
%     t_vec = 0:dt:out.t(end); 
% end
% 
% [tout,Dout] = ode45(@(t,D) dataTime(t,D,phase_type,debrisRateCoeff,out.t,lat,DLIntervals,dataPerDebris,K_DL,R_GNC,R_EPS,R_Struct,pi_cVal), ...
%     t_vec,[0; 0],opts);
% 
% debrisRateOut = interp1(out.t,debrisRate,tout)'; % events/second

t_vec = 0:dt:out.t(end);
nSteps = numel(t_vec);
D = zeros(2, nSteps);
dDdt = zeros(1, nSteps);

for k = 1:nSteps
    
    f = dataTime(t_vec(k),D(:,k),phase_type,debrisRateCoeff,out.t,lat,DLIntervals,dataPerDebris,K_DL,R_GNC,R_EPS,R_Struct,pi_cVal);
    dDdt(k) = f(1);
    if k < nSteps
        D(:,k+1) = max(D(:,k) + dt*f, 0);   % stored data can't go negative
    end
end
tout= t_vec';
Dout = D';

% DLIdx = any(cell2mat(arrayfun(@(i) tout>=DLIntervals(1,i) & tout<=DLIntervals(2,i),1:length(DLIntervals),'UniformOutput',false))');
% debrisRates = (debrisRateCoeff*normpdf(interp1(out.t,lat,tout),0,90))'; 
% R_payload = debrisRates * dataPerDebris .* double(~DLIdx);
% 
% % only downlink data that you can downlink
% K_DLmin = min(Dout(:,1) + 1/pi_cVal*(R_GNC+R_EPS+R_Struct),K_DL)' .* double(DLIdx); 
% 
% dDdt = 1/pi_cVal*(R_payload+R_GNC+R_EPS+R_Struct) - K_DLmin;

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

function dDdt = dataTime(t,D, phase_type, debrisRateCoeff,oldTime,lat,DLIntervals,dataPerDebris,K_DL,R_GNC,R_EPS,R_Struct,pi_c)


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
