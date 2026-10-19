addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..'));
odinPath;
res = tradeLighting(365);
plotTradeLighting(res, 5);