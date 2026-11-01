clear; clc; close all;

rawData = readmatrix("synthetic_debris_events.txt");
N = length(rawData);

% posX = uint64(rawData(:,1));
% posY = uint64(rawData(:,1));
% pol = uint64(rawData(:,1));
% timestamp = uint64(rawData(:,1));
% 
% results = bitshift(posX, 52) + ...
%              bitshift(posY, 40) + ...
%              bitshift(pol,  39) + ...
%              timestamp;

% 1. Force everything to uint64 immediately
posX      = uint64(rawData(:,1)); 
posY      = uint64(rawData(:,2)); 
pol       = uint64(rawData(:,3)); 
timestamp = uint64(rawData(:,4)); 

% 2. Cleanly mask each value to its exact bit-width to clear any ghost bits
posX_masked      = bitand(posX,      uint64(2^12 - 1)); % Keep lower 12 bits
posY_masked      = bitand(posY,      uint64(2^12 - 1)); % Keep lower 12 bits
pol_masked       = bitand(pol,       uint64(2^1 - 1));  % Keep 1 bit
timestamp_masked = bitand(timestamp, uint64(2^39 - 1)); % Keep lower 39 bits

% 3. Combine safely using bitor instead of addition (+)
packedData = bitor(bitshift(posX_masked, 52), ...
             bitor(bitshift(posY_masked, 40), ...
             bitor(bitshift(pol_masked,  39), ...
             timestamp_masked)));
%% write to file
f = fopen("synthetic_debris_events_noise.bin","wb");
fwrite(f,packedData,'uint64');
fclose(f);

f = fopen("synthetic_debris_events_noise.bin","r");
encData = fread(f,Inf,"uint64=>uint64");
fclose(f);

assert(all(encData==packedData),"Data does not match!");