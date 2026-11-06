function frameDrift(nYear)
% frameDrift: validation of the model frame over the mission, the Earth fixed directions from earthRotation.m against the Aerospace Toolbox IAU-2000/2006 Earth frame of date.
% nYear: run length, 5 when omitted [years]

if nargin < 1
    nYear = 5;
end

here = fileparts(mfilename('fullpath')); % here: folder holding this file
root = fullfile(here, '..'); % root: repository root folder
addpath(root);
odinPath;

outDir = fullfile(here, 'results'); % outDir: folder for logs and results
if ~isfolder(outDir)
    mkdir(outDir);
end
logFile = fullfile(outDir, sprintf('frameDrift_%dy.txt', nYear)); % logFile: text log of this run
if isfile(logFile)
    delete(logFile);
end
diary(logFile);

cfg = config(); % cfg: parameter struct from config.m
nPt = 200*nYear; % nPt: number of random points, 200 a year
lim = 0.01; % lim: limit on the largest separation over the run, as declared for the astropy check in Checks 6 [deg]
dAT = 37; % dAT: TAI minus UTC, as in gcrfRotation.m [s]
yr = 365.25*86400; % yr: one Julian year [s]

rng(1, 'twister');
t = sort(nYear*yr*rand(1, nPt)); % t: 1 x nPt times since sim start, uniform over the run [s]
u = randn(3, nPt); % u: 3 x nPt directions in the model frame, uniform over the sphere
u = u./vecnorm(u);

th = earthRotation(t, cfg); % th: 1 x nPt Greenwich angle [rad]

% Model frame to Earth fixed, rotation about z by the Greenwich angle, as in dragDensity.m:
eM = [cos(th).*u(1,:) + sin(th).*u(2,:); cos(th).*u(2,:) - sin(th).*u(1,:); u(3,:)]; % eM: 3 x nPt Earth fixed directions from the model

R = gcrfRotation(cfg); % R: model frame to GCRF rotation
utc = datevec(datenum(cfg.epoch) + t(:)/86400); % utc: nPt x 6 UTC date vectors
D = dcmeci2ecef('IAU-2000/2006', utc, dAT*ones(nPt, 1)); % D: 3 x 3 x nPt GCRF to Earth fixed rotations, UT1 = UTC, no polar motion

eT = zeros(3, nPt); % eT: 3 x nPt Earth fixed directions from the Toolbox
for k = 1:nPt
    eT(:,k) = D(:,:,k)*R*u(:,k);
end

% Angle between the two Earth fixed directions:
sep = atan2d(vecnorm(cross(eM, eT)), dot(eM, eT)); % sep: 1 x nPt ground point separation [deg]

ty = (0:nYear)*yr; % ty: 1 x (nYear + 1) start of each year and end of the run [s]
Dy = dcmeci2ecef('IAU-2000/2006', datevec(datenum(cfg.epoch) + ty(:)/86400), dAT*ones(nYear + 1, 1)); % Dy: 3 x 3 x (nYear + 1) GCRF to Earth fixed rotations
pole = zeros(1, nYear + 1); % pole: 1 x (nYear + 1) angle of the model pole from the pole of date [deg]
for k = 1:nYear + 1
    z = Dy(:,:,k)*R*[0; 0; 1]; % z: model pole in the Earth fixed frame of date
    pole(k) = atan2d(hypot(z(1), z(2)), z(3));
end

fprintf('model frame, mean equator and equinox of %04d-%02d-%02d %02d:%02d UTC through earthRotation.m, against gcrfRotation.m then dcmeci2ecef IAU-2000/2006\n', ...
        cfg.epoch(1:5));
fprintf('UT1 = UTC, no polar motion, TAI - UTC %d s; %d random directions, uniform over %d years, seed 1\n\n', dAT, nPt, nYear);
fprintf('year   points   largest ground point separation [deg]   pole of epoch from the pole of date at the start [deg]\n');
for y = 1:nYear
    in = t >= (y - 1)*yr & t < y*yr; % in: points in year y
    fprintf('%4d   %6d   %37.4f   %54.4f\n', y, nnz(in), max(sep(in)), pole(y));
end
fprintf('end of year %d: pole of epoch from the pole of date %.4f deg\n', nYear, pole(end));
ok = max(sep) <= lim; % ok: pass flag of the separation over the run
fprintf('largest separation over %d years %.4f deg, limit %.2f deg: %s\n', nYear, max(sep), lim, verdict(ok));

diary off;
save(fullfile(outDir, sprintf('frameDrift_%dy.mat', nYear)), 'sep', 't', 'u', 'pole', 'ok');

end


function s = verdict(pass)
% verdict: PASS or CHECK for a logical test result.
% pass: true when the check passed
% s: text for the log

if pass
    s = 'PASS';
else
    s = 'CHECK';
end

end
