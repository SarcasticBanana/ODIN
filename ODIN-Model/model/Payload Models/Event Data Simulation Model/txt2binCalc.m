%% Custom 64-bit word 
clear; clc; close all;

rawData = readmatrix("synthetic_debris_events.txt");
N = length(rawData);
nBits = N*64;
tTot = rawData(end,4)/1e6;

bps = nBits/tTot;
disp(['Output: ', num2str(bps/1e6), ' Mbps']);

%% EVT 3.0 16-bit word encoding
clear; clc; close all;

rawData = readmatrix("synthetic_debris_events.txt");
N = length(rawData);
xs = rawData(:,1);
ys = rawData(:,2);
pols = rawData(:,3);
times = rawData(:,4);
clear rawData;

%%
nBitsTime = 32*numel(unique(times));

reprepYs = numel(find(ys==circshift(ys,1) & ys==circshift(ys,2)));
repYs = numel(find(ys==circshift(ys,1))) - reprepYs*2;
uniqueYs = numel(ys)-repYs;
nBitsY = 16*(uniqueYs - repYs - 2*reprepYs);
nBitsX = 16*(numel(xs));

nBitsTot = nBitsTime + nBitsY + nBitsX;
bps = nBitsTot/(times(end)/1e6);
disp(['Output: ', num2str(bps/1e6), ' Mbps']);
