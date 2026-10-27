%%
clear; clc; close all;

N = 16*1e4;
headerLen = 316*8;
file = fopen("DUBHE_2025-02-26_19-37-34.raw");
rawDataBytes = fread(file,N+headerLen,"uint8");
fclose(file);
bits_matrix = de2bi(rawDataBytes, 8, 'right-msb');
rawData = reshape(bits_matrix.', [], 1);

% cut off header
rawData = rawData(headerLen+1:end);

% group into 16-bit "event types"
nBits = 16;
groupData = fliplr(reshape(rawData,nBits,[])');
strGroups = arrayfun(@(i) strjoin(string(groupData(i,:)),""),1:N/nBits)';

% name event types
strGroups(startsWith(strGroups,"0000"),2) = "EVT_ADDR_Y";
strGroups(startsWith(strGroups,"0010"),2) = "EVT_ADDR_X";
strGroups(startsWith(strGroups,"0011"),2) = "VECT_BASE_X";
strGroups(startsWith(strGroups,"0100"),2) = "VECT_12";
strGroups(startsWith(strGroups,"0101"),2) = "VECT_8";
strGroups(startsWith(strGroups,"0110"),2) = "EVT_TIME_LOW";
strGroups(startsWith(strGroups,"0111"),2) = "CONTINUED_4";
strGroups(startsWith(strGroups,"1000"),2) = "EVT_TIME_HIGH";
strGroups(startsWith(strGroups,"1010"),2) = "EXT_TRIGGER";
strGroups(startsWith(strGroups,"1110"),2) = "OTHERS";
strGroups(startsWith(strGroups,"1111"),2) = "CONTINUED_12";
%%
EVT_TIME_LOW = char(strGroups(startsWith(strGroups,"0110")));
EVT_TIME_LOW = EVT_TIME_LOW(:,5:end);
EVT_TIME_HIGH = char(strGroups(startsWith(strGroups,"1000")));
EVT_TIME_HIGH = EVT_TIME_HIGH(:,5:end);
% EVT_TIME_HIGH should be stored in variable until it updates

EVT_TIME_HIGH = char(strGroups(find(startsWith(strGroups,"1000"),1)));
EVT_TIME_HIGH = EVT_TIME_HIGH(5:end);
EVT_TIME_LOW = char(strGroups(find(startsWith(strGroups,"0110"),1)));
EVT_TIME_LOW = EVT_TIME_LOW(5:end);

EVT_TIME(1) = bin2dec([EVT_TIME_HIGH,EVT_TIME_LOW]);
% dataFrame = zeros(1280,720,length(strGroups));
dataFrame = zeros(1280,720);
frame = zeros(1280,720);

figure;
plt = imagesc(dataFrame(:,:,1)');

colorbar;
axis image;
for i=2:length(strGroups)
    frame = zeros(1280,720);
    switch extractBefore(strGroups(i),5)
        case "1000"
            EVT_TIME_HIGH = char(strGroups(i));
            EVT_TIME_HIGH = EVT_TIME_HIGH(5:end);
        case "0110"
            EVT_TIME_LOW = char(strGroups(i));
            EVT_TIME_LOW = EVT_TIME_LOW(5:end);
        case "0000"
            posYdata = char(strGroups(i));
            posY = bin2dec(posYdata(6:end));
        case "0010"
            posXdata = char(strGroups(i));
            posX = bin2dec(posXdata(6:end));
            polarity = bin2dec(posXdata(5))*2 - 1;
            frame(posX+1,posY+1) = polarity;
        otherwise
    end

    EVT_TIME(i) = bin2dec([EVT_TIME_HIGH,EVT_TIME_LOW]);
    if any(frame,'all')
        % dataFrame(:,:,i) = dataFrame(:,:,i-1) + frame;
        dataFrame = dataFrame + frame;
        % disp(['changing pixel at ',num2str(posX+1),',',num2str(posY+1),' in ',num2str(polarity),' direction']);
    % else
    %     dataFrame(:,:,i) = dataFrame(:,:,i-1);
    end

    % plt.CData = dataFrame(:,:,i)';
    if mod(i,250)==0
        plt.CData = dataFrame';
        % xlim([900 1000]);
        % ylim([400 500]);
        title(['t = ',num2str(EVT_TIME(i)/1e6)]);
        drawnow limitrate;
    end
end
