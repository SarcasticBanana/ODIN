clear; close all; clc;

addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
odinPath;

cfg = config();
cfg.h = 600;                  % Your chosen altitude [km]
cfg.ode.dtOut = 1*60;
nOrbit = 150; % change this for 10 day sim
out = propagate(cfg, nOrbit);     % Ale's orbit calculation
N = length(out.t);

phase_type = "Calibrate"; % Phase 1 = "Calibrate", Phase 2 = "Nominal"

% lat = sign(out.r(3,:)) .* acosd(vecnorm(out.r(1:2,:)) ./ vecnorm(out.r));
% figure; plot(out.t,lat);
% figure; plot(lat,normpdf(lat,0,120));

% [rate_per_s,details] = debrisDetectionRate(out.t,out.r',out.v',vecnorm(out.r)'-cfg.Re,cfg,out.rsun');
% debrisRate = ones(1,N)*details.TotalRate_per_s;

load('debrisRates.mat');
debrisT = xData;
debrisRates = yData;
debrisRates = debrisRates(out.t<=out.T);
debrisRate = repmat(debrisRates,1,nOrbit);

dataPerDebris = 22e6;
% dataBG = 17e6;
% dataDebris = 20e6; % when detection

% t = 3600*24*10;

st = stations();
[dl,c] = downlinkTimes(out,cfg,st);
DLIntervals = dl';

K_DL = 56e6; % 56 Mbps, probably
R_GNC = 100; % 100 b/s, guess
R_EPS = 100;
R_Struct = 10; % 10 b/s, guess
headerBits = 48;
dataBitsPerHeader = 65536*8;
packetizationInflation = (headerBits+dataBitsPerHeader)/dataBitsPerHeader; 
pi_c = 2.34/packetizationInflation; % compression ratio, probably. Note: should be division not mulitplication? the inflation should decrease amount of compression?
storage_cap = 64 * 8e9; % [b] Upper bound on storage SD card
recollect_threshold = 0; % [b] Describes threshold required during calibration to start re-collecting data after SD card fills up.

dt = 0.5;
global eventHappening tEvent
eventHappening = 0;
tEvent = 0;

if strcmpi(phase_type, "Calibrate")
    opts = odeset('NonNegative',1,'MaxStep',10); % just collect data, no need to model probabilist events
    t_vec = [0, out.t(end)]; % Let ode45 choose step
else
    opts = odeset('NonNegative',1,'MaxStep',10,'OutputFcn', @odeOutput);
    t_vec = 0:dt:out.t(end); 
end

[tout,Dout] = ode45(@(t,D) dataTime(t,D,phase_type,DLIntervals,dataPerDebris,K_DL,R_GNC,R_EPS,R_Struct,pi_c, storage_cap), ...
    t_vec,[0; 0],opts);
%
% debrisRateOut = interp1(out.t,debrisRate,tout)'; % events/second

DLIdx = any(cell2mat(arrayfun(@(i) tout>=DLIntervals(1,i) & tout<=DLIntervals(2,i),1:length(DLIntervals),'UniformOutput',false))');

R_payload =  (dataPerDebris) .* double(~DLIdx);

% only downlink data that you can downlink
K_DLmin = min(Dout(:,1) + 1/pi_c*(R_GNC+R_EPS+R_Struct),K_DL)' .* double(DLIdx); 

dDdt = 1/pi_c*(R_payload+R_GNC+R_EPS+R_Struct) - K_DLmin;
%%
figure; hold on;
plot(tout,Dout(:,1),'b','LineWidth',2);
ylabel("Data Stored (Bits)");
yline(storage_cap, '-b')
yyaxis('right');
% plot(tout,dDdt,'r','LineWidth',2);
yline(0,'r--');
% xline(out.T*(1:nOrbit),'k--');
xline(24*3600*(1:out.t(end)/(24*3600)));
ylabel("Total Data Rate (Bits/s)");
xlabel("Time (s)");
box on; grid on;
ax = gca;
ax.YAxis(1).Color = [0 0 1];
ax.YAxis(2).Color = [1 0 0];

function status = odeOutput(t,D,flag)
    global eventHappening tEvent
    if isempty(flag)
        if ~eventHappening && length(t)>1
            % debrisRate = interp1(oldTime,debrisRate,t(end));
            debrisRate = 1.98e-3; % assume const avg debris rate
            dt = t(end)-t(1);
            % pEvent = exp(-debrisRate)*debrisRate^(dt);
            pEvent = debrisRate*dt*exp(-debrisRate*dt);
            eventHappening = rand(1)<pEvent;
            if eventHappening
                % disp(['probability of event happening in ',num2str(dt),' is ',num2str(pEvent)]);
                tEvent = t(end);
                disp(['time %: ',num2str(t(end)/290040*100)]);
            end
        elseif eventHappening && (t(end)-tEvent)>randn(1)*5+8
            eventHappening = 0;
        end
        
        status = 0;
    end
end

function dDdt = dataTime(t,D,phase_type,DLIntervals,dataPerDebris,K_DL,R_GNC,R_EPS,R_Struct,pi_c, storage_cap)
    global eventHappening
    DLIdx = any(arrayfun(@(i) t>=DLIntervals(1,i) & t<=DLIntervals(2,i),1:length(DLIntervals)));

    if strcmpi(phase_type,"Calibrate")
        R_payload = dataPerDebris*double(~DLIdx); % for calibration
    elseif strcmpi(phase_type,"Nominal")
        R_payload = dataPerDebris * double(~DLIdx) * double(eventHappening);
    end
    % TODO: payload should also be stopped in eclipse
    if D(1) >= storage_cap
        R_payload = 0; % Storage is full, camera is stopped.
    end

    K_DLmin = min(D(1) + 1/pi_c*(R_GNC + R_EPS + R_Struct),K_DL) * double(DLIdx); % only downlink data that you can downlink for min bound

    dDdt = [1/pi_c*(R_payload+R_GNC+R_EPS+R_Struct) - K_DLmin;
            0];
end
