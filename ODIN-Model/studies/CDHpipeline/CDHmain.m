clear; close all; clc;

addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
odinPath;

cfg = config();
cfg.h = 600;                  % Your chosen altitude [km]
cfg.ode.dtOut = 1*60;
nOrbit = 15; % change this for 10 day sim
out = propagate(cfg, nOrbit);     % Ale's orbit calculation
N = length(out.t);

phase_type = "Calibrate"; % Phase 1 = "Calibrate", Phase 2 = "Nominal"

% [rate_per_s,details] = debrisDetectionRate(out.t,out.r',out.v',vecnorm(out.r)'-cfg.Re,cfg,out.rsun');
% debrisRate = ones(1,N)*details.TotalRate_per_s;

load('debrisRates.mat');
debrisT = xData;
debrisRates = yData;
debrisRates = debrisRates(out.t<=out.T);
debrisRate = repmat(debrisRates,1,nOrbit);

noiseSamples = load('model/Payload Models/Event Data Simulation Model/eventSamplesNoise.mat');
dataRateTimes = noiseSamples.tspost;
dataRates = noiseSamples.dataRate64;

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
storage_cap = 0.8*64 * 8e9; % [b] Upper bound on storage SD card
recollect_threshold = 0; % [b] Describes threshold required during calibration to start re-collecting data after SD card fills up.

dt = 0.5;
global eventHappening tEvent downlinking
eventHappening = 0;
tEvent = 0;
downlinking = 0;

if strcmpi(phase_type, "Calibrate")
    opts = odeset('NonNegative',1,'MaxStep',10,'AbsTol',1,'OutputFcn', @odeOutput); % just collect data, no need to model probabilist events
    t_vec = [0, out.t(end)]; % Let ode45 choose step
else
    opts = odeset('NonNegative',1,'MaxStep',10,'AbsTol',1,'OutputFcn', @odeOutput);
    t_vec = 0:dt:out.t(end); 
end

[tout,Dout] = ode45(@(t,D) dataTime(t,D,phase_type,dataRateTimes,dataRates,DLIntervals,K_DL,R_GNC,R_EPS,R_Struct,pi_c,storage_cap), ...
    t_vec,[0; 0],opts);
%
% debrisRateOut = interp1(out.t,debrisRate,tout)'; % events/second

DLIdx = any(cell2mat(arrayfun(@(i) tout>=DLIntervals(1,i) & tout<=DLIntervals(2,i),1:length(DLIntervals),'UniformOutput',false))');

% R_payload =  dataRates .* double(~DLIdx);
% 
% % only downlink data that you can downlink
K_DLmin = min(Dout(:,1) + 1/pi_cAvg*(R_GNC+R_EPS+R_Struct),K_DL)' .* double(DLIdx); 
% 
% dDdt = 1/pi_c*(R_payload+R_GNC+R_EPS+R_Struct) - K_DLmin;

disp(['Max Data Stored: ',num2str(max(Dout(:,1))/1e9),'Gb = ',num2str(max(Dout(:,1))/8e9),'GB of data!']);
totDownlinked = trapz(tout,K_DLmin)*pi_cAvg;
disp(['Downlinked ',num2str(totDownlinked/1e9),'Gb = ',num2str(totDownlinked/8e9),'GB of data!']);

figure; hold on;
plot(tout/24/3600,Dout(:,1)/1e9);
ylabel("Data Stored (Gb)");
% yline(16*8-0.5*8)
% yline(storage_cap/1e9, '-b');
% yyaxis('right');
% plot(tout,dDdt,'r','LineWidth',2);
% yline(0,'r--');
% xline(out.T*(1:nOrbit),'k--');
% xline(24*3600*(1:out.t(end)/(24*3600)));
xlim([0 out.t(end)/24/3600])
% ylabel("Total Data Rate (Bits/s)");
xlabel("Time (days)");
box on; grid on;
ax = gca;
% ax.YAxis(1).Color = [0 0 1];
% ax.YAxis(2).Color = [1 0 0];

function status = odeOutput(t,D,flag)
    global eventHappening tEvent
    if isempty(flag)
        % if ~eventHappening && length(t)>1
        %     % debrisRate = interp1(oldTime,debrisRate,t(end));
        %     debrisRate = 1.98e-3; % assume const avg debris rate
        %     dt = t(end)-t(1);
        %     % pEvent = exp(-debrisRate)*debrisRate^(dt);
        %     pEvent = debrisRate*dt*exp(-debrisRate*dt);
        %     eventHappening = rand(1)<pEvent;
        %     if eventHappening
        %         % disp(['probability of event happening in ',num2str(dt),' is ',num2str(pEvent)]);
        %         tEvent = t(end);
        %         disp(['time %: ',num2str(t(end)/870180*100)]);
        %     end
        % elseif eventHappening && (t(end)-tEvent)>randn(1)*5+8
        %     eventHappening = 0;
        % end
        if any(isapprox(t(end)/870180,0.001:0.001:1,'loose'))
            disp(['time %: ',num2str(t(end)/870180*100)]);
        end
        status = 0;
    end
end

function dDdt = dataTime(t,D,phase_type,dataRateTimes,dataRates,DLIntervals,K_DL,R_GNC,R_EPS,R_Struct,pi_c,storage_cap)
    global eventHappening

    DLIdx = any(arrayfun(@(i) t>=DLIntervals(1,i) & t<=DLIntervals(2,i),1:length(DLIntervals)));
    ObservingIdx = ~any(arrayfun(@(i) t>=DLIntervals(1,i)-300 & t<=DLIntervals(2,i)+300,1:length(DLIntervals)));

    % if D(1)>=storage_cap
    %     downlinking = 1;
    % elseif downlinking && D(1)<10
    %     downlinking = 0;
    % end

    if strcmpi(phase_type,"Calibrate")
        if t>max(dataRateTimes)
            dataRate = 15e6;
        else
            dataRate = interp1(dataRateTimes,dataRates,t,'linear','extrap');
        end
        R_payload = dataRate*double(ObservingIdx); % for calibration
    elseif strcmpi(phase_type,"Nominal")
        if t>max(dataRateTimes)
            dataRate = 15e6;
        else
            dataRate = interp1(dataRateTimes,dataRates,t,'linear','extrap');
        end
        R_payload = dataRate * double(ObservingIdx) * double(eventHappening);
    end

    % if D(1) >= storage_cap
    %     R_payload = 0; % Storage is full, camera is stopped.
    % end

    thisPi_c = pi_c(1);
    K_DLmin = min(D(1) + 1/thisPi_c*(R_GNC + R_EPS + R_Struct),K_DL) * double(DLIdx); % only downlink data that you can downlink for min bound

    dDdt = [1/thisPi_c*(R_payload+R_GNC+R_EPS+R_Struct) - K_DLmin;
            0];
end
