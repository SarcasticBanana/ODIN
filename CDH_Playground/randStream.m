clear; clc; close all;

n = 1e3;

posX = round(rand(n,11));
posY = round(rand(n,11));
bit = round(rand(n,1));
time = int2bit(1:n,32)';

output = arrayfun(@(i) strjoin(string([posX(i,:),posY(i,:),bit(i),time(i,:)]),""),1:n)';
%%
clear; clc; close all;

N = 16e3;
headerLen = 316*8;
file = fopen("DUBHE_2025-02-26_19-37-34.raw");
rawData = fread(file,N+headerLen,"ubit1");
fclose(file);
rawData = rawData(headerLen+1:end);

nBits = 16;
groupData = fliplr(reshape(rawData,nBits,[])');
strGroups = arrayfun(@(i) strjoin(string(groupData(i,:)),""),1:N/nBits)';

EVT_TIME_LOW = char(strGroups(startsWith(strGroups,"0110")));
EVT_TIME_LOW = EVT_TIME_LOW(:,5:end);
EVT_TIME_HIGH = char(strGroups(startsWith(strGroups,"1000")));
EVT_TIME_HIGH = EVT_TIME_HIGH(:,5:end);

%%
figure;
% histogram(rawData,"BinWidth",1);
findEvents = find(rawData==1 | rawData==0);
distEvents = findEvents-circshift(findEvents,1);
distEvents(1) = nan;
histogram(distEvents,"BinWidth",1);
% xlim([0 32]);
xlabel("Distance between 1 or 0");