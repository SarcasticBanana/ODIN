clear; clc; close all;
odinPath();

nx = 2560;
ny = 2560;
N = 20;
% odd indices are x=0, y=any (with pos x vel)
% even indices are x=any, y=0 (with pos y vel)
x0s = round(rand(1,N)*nx);
x0s(1:2:end) = 0;
y0s = round(rand(1,N)*ny);
y0s(2:2:end) = 0;
vMax = 3e3; % for a max diagonal event, the fastest an object can go to cross in 0.5s
vxs = randn(1,N)*vMax;
vxs(2:2:end) = vxs(2:2:end).*(round(rand(1,N/2))*2-1);
vys = rand(1,N)*vMax;
vys(1:2:end) = vys(1:2:end).*(round(rand(1,N/2))*2-1);

ts = min(abs((nx.*(vxs>0)-x0s)./vxs),abs((ny.*(vys>0)-y0s)./vys));
tfilter = ts<0.5 | ts>20;
x0s(tfilter) = [];
y0s(tfilter) = [];
vxs(tfilter) = [];
vys(tfilter) = [];
N = nnz(~(tfilter));
disp(['Proceeding with ', num2str(N), ' long-enough events']);

figure; quiver(x0s,y0s,vxs,vys);
xline([0,nx]);
yline([0,nx]);
%% using random x0, y0 values
% implement phi_0 photon flux 
Phi0 = 688000;
for i=1:N
    [~] = eventDataStreamFun(x0s(i),y0s(i),vxs(i),vys(i),Phi0);
end

%% full random noise
clear; clc; close all;
odinPath();

nx = 2560;
ny = 2560;
N = 1;

Ts = linspace(1,30,N);

Phi0 = 0;

for i=1:N
    [~] = eventDataStreamFun(-nx,-ny,0,0,Phi0,Ts(i));
end

%%
% filelist = dir("model/Payload Models/Event Data Simulation Model/census/*.txt");
filelist = dir("model/Payload Models/Event Data Simulation Model/censusNoise/*.txt");
oldSamples = load("eventSamplesNoise.mat");
Npost = length(filelist);
dataRate64 = zeros(1,Npost);
dataRate16 = zeros(1,Npost);
detectable = zeros(1,Npost);
x0post = zeros(1,Npost);
y0post = zeros(1,Npost);
vxpost = zeros(1,Npost);
vypost = zeros(1,Npost);
tspost = zeros(1,Npost);

for i=1:Npost
    filepath = [filelist(i).folder,'\',filelist(i).name];
    [dataRate64(i),dataRate16(i)] = txt2binCalc(filepath);
    if ~contains(filelist(i).folder,'Noise')
        detectable(i) = AngularRateFun(filepath);
    end

    regex = regexp(filelist(i).name,"_x(.+)_y(.+)_vx(.+)_vy(.+)_tsim(.+).txt",'tokens');
    x0post(i) = str2double(regex{1}{1});
    y0post(i) = str2double(regex{1}{2});
    vxpost(i) = str2double(regex{1}{3});
    vypost(i) = str2double(regex{1}{4});
    tspost(i) = str2double(regex{1}{5});
end

Npost = Npost+oldSamples.Npost;
dataRate64 = [oldSamples.dataRate64,dataRate64];
dataRate16 = [oldSamples.dataRate16,dataRate16];
detectable = [oldSamples.detectable,detectable];
x0post = [oldSamples.x0post,x0post];
y0post = [oldSamples.y0post,y0post];
vxpost = [oldSamples.vxpost,vxpost];
vypost = [oldSamples.vypost,vypost];
tspost = [oldSamples.tspost,tspost];
% tspost = min(abs((nx.*(vxpost>0)-x0post)./vxpost),abs((ny.*(vypost>0)-y0post)./vypost));

% save model\'Payload Models'\'Event Data Simulation Model'\eventSamples.mat Npost x0post y0post vxpost vypost dataRate detectable tspost
save model\'Payload Models'\'Event Data Simulation Model'\eventSamplesNoise.mat Npost x0post y0post vxpost vypost dataRate64 dataRate16 detectable tspost

%%
close all;
% load model\'Payload Models'\'Event Data Simulation Model'\eventSamples.mat
% load model\'Payload Models'\'Event Data Simulation Model'\eventSamplesNoise.mat

figure;
quiver(x0post,y0post,vxpost,vypost);
grid on;
xline([0,nx])
yline([0,nx])

figure;
histogram(vecnorm([vxpost;vypost]),round(Npost/4));
xlabel("Magnitude of Velocity (px/s)");

figure;
histogram(dataRate64/1e6,round(Npost/4));
xlabel('Data Rate (Mbps)');

detectLog = logical(detectable);
figure; hold on;
scatter(vecnorm([vxpost(detectLog);vypost(detectLog)]),dataRate64(detectLog)/1e6,'b');
scatter(vecnorm([vxpost(~detectLog);vypost(~detectLog)]),dataRate64(~detectLog)/1e6,'r');
grid on;
legend("Detectable","Not Detectable");
xlabel('Magnitude of Velocity (px/s)');
ylabel('Data Rate (Mbps)');

figure; hold on;
if all(~detectLog)
    scatter(tspost,dataRate64/1e6);
else
    scatter(tspost(detectLog),dataRate64(detectLog)/1e6,'b');
    scatter(tspost(~detectLog),dataRate64(~detectLog)/1e6,'r');
    legend("Detectable","Not Detectable");
end
grid on;
xlabel('Time in FOV (s)');
ylabel('Data Rate (Mbps)');

%%
function bps = fitDataRate(filename)
    rawData = readmatrix(filename);
    N = length(rawData);
    times = rawData(:,4);
    % No10 = round(N/100);
    % [p,S] = polyfit(times(N-No10:N)/1e6,0:No10,1);
    % [fitresult,gof] = fit(times(N-No10:N),(0:No10)','poly1');
    % eventsPerSec = fitresult.p1*1e6;
    % bps = eventsPerSec*64;
    % fprintf("Fit to slope of %.0f events/sec = %.2f Mbps with R^2 = %.6f\n",eventsPerSec,bps/1e6,gof.rsquare);

    binW = 0.25e6;
    centerBins = binW/2:binW:times(end);
    tbins = [centerBins-binW/2;
             centerBins+binW/2];
    for i=1:length(centerBins)
        tInBin = times(times>tbins(1,i) & times<tbins(2,i));
        bpsAtT(i) = numel(tInBin)*1e6/(tInBin(end)-tInBin(1));
    end
    % figure; plot(centerBins,bpsAtT)
    % xlabel('Time');
    % ylabel('Data Rate (bps)');
    % figure; plot(1./sqrt(centerBins),bpsAtT)
    % xlabel('Time (1/sqrt(s))');
    % ylabel('Data Rate (bps)');
    ft = fittype( 'a*x*exp(c*x)+b*exp(d*x)+e', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Lower = [0 -Inf -Inf -Inf 0];
    opts.StartPoint = [0.05, -0.5, -0.5, -0.1, 0.5];
    opts.Upper = [Inf 0 0 0 Inf];
    [fitresult, gof] = fit((centerBins/1e6)', (bpsAtT/1e6)', ft, opts);
    eventsPerSec = fitresult.e*1e6;
    bps = eventsPerSec*64;

    figure;
    plot(fitresult, (centerBins/1e6)', (bpsAtT/1e6)');
    grid on;


    fprintf("Fit to final slope of %.0f events/sec = %.2f Mbps with R^2 = %.6f\n",eventsPerSec,bps/1e6,gof.rsquare);
end

function [bps64,bps16] = txt2binCalc(filename)
    % 64-bit encoding
    rawData = readmatrix(filename);
    N = length(rawData);
    nBits = N*64;
    tTot = rawData(end,4)/1e6;
    
    bps64 = nBits/tTot;
    disp(['64-bit Output: ', num2str(bps64/1e6), ' Mbps']);

    % 16-bit encoding
    xs = rawData(:,1);
    ys = rawData(:,2);
    times = rawData(:,4);
    clear rawData;
    
    nBitsTime = 32*numel(unique(times));
    
    reprepYs = numel(find(ys==circshift(ys,1) & ys==circshift(ys,2)));
    repYs = numel(find(ys==circshift(ys,1))) - reprepYs*2;
    uniqueYs = numel(ys)-repYs;
    nBitsY = 16*(uniqueYs - repYs - 2*reprepYs);
    nBitsX = 16*(numel(xs));
    
    nBitsTot = nBitsTime + nBitsY + nBitsX;
    bps16 = nBitsTot/(times(end)/1e6);
    disp(['16-bit Output: ', num2str(bps16/1e6), ' Mbps']);
end