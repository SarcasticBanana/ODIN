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
vMax = 1e3; % for a max diagonal event, the fastest an object can go to cross in 0.5s
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
%%
for i=1:N
    [events] = eventDataStreamFun(x0s(i),y0s(i),vxs(i),vys(i));
end
%%
filelist = dir("model/Payload Models/Event Data Simulation Model/census/*.txt");
Npost = length(filelist);
dataRate = zeros(1,Npost);
detectable = zeros(1,Npost);
x0post = zeros(1,Npost);
y0post = zeros(1,Npost);
vxpost = zeros(1,Npost);
vypost = zeros(1,Npost);

for i=1:Npost
    filepath = [filelist(i).folder,'\',filelist(i).name];
    dataRate(i) = txt2binCalc(filepath);
    detectable(i) = AngularRateFun(filepath);

    regex = regexp(filelist(i).name,"_x(.+)_y(.+)_vx(.+)_vy(.+).txt",'tokens');
    x0post(i) = str2double(regex{1}{1});
    y0post(i) = str2double(regex{1}{2});
    vxpost(i) = str2double(regex{1}{3});
    vypost(i) = str2double(regex{1}{4});
end
tspost = min(abs((nx.*(vxpost>0)-x0post)./vxpost),abs((ny.*(vypost>0)-y0post)./vypost));

save model\'Payload Models'\'Event Data Simulation Model'\eventSamples.mat Npost x0post y0post vxpost vypost dataRate detectable tspost
%%
figure;
quiver(x0post,y0post,vxpost,vypost);
grid on;
xline([0,nx])
yline([0,nx])

figure;
histogram(vecnorm([vxpost;vypost]),round(Npost/4));
xlabel("Magnitude of Velocity (px/s)");

figure;
histogram(dataRate/1e6,round(Npost/4));
xlabel('Data Rate (Mbps)');

detectLog = logical(detectable);
figure; hold on;
scatter(vecnorm([vxpost(detectLog);vypost(detectLog)]),dataRate(detectLog)/1e6,'b');
scatter(vecnorm([vxpost(~detectLog);vypost(~detectLog)]),dataRate(~detectLog)/1e6,'r');
grid on;
legend("Detectable","Not Detectable");
xlabel('Magnitude of Velocity (px/s)');
ylabel('Data Rate (Mbps)');

figure; hold on;
scatter(tspost(detectLog),dataRate(detectLog)/1e6,'b');
scatter(tspost(~detectLog),dataRate(~detectLog)/1e6,'r');
grid on;
legend("Detectable","Not Detectable");
xlabel('Time in FOV (s)');
ylabel('Data Rate (Mbps)');


%% EVT 3.0 16-bit word encoding
function bps = txt2binCalc(filename)
    rawData = readmatrix(filename);
    % N = length(rawData);
    xs = rawData(:,1);
    ys = rawData(:,2);
    % pols = rawData(:,3);
    times = rawData(:,4);
    clear rawData;
    
    nBitsTime = 32*numel(unique(times));
    
    reprepYs = numel(find(ys==circshift(ys,1) & ys==circshift(ys,2)));
    repYs = numel(find(ys==circshift(ys,1))) - reprepYs*2;
    uniqueYs = numel(ys)-repYs;
    nBitsY = 16*(uniqueYs - repYs - 2*reprepYs);
    nBitsX = 16*(numel(xs));
    
    nBitsTot = nBitsTime + nBitsY + nBitsX;
    bps = nBitsTot/(times(end)/1e6);
    disp(['Output: ', num2str(bps/1e6), ' Mbps']);
end