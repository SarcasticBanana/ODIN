%% EVENT_DATA_DIAGNOSTICS
% Looks at an event file (x  y  polarity  t_us) and prints / plots the numbers
% you need to calibrate the simulator against real data.
%
% Run this SAME script on:
%   (1) the real event file   -> these are your calibration TARGETS
%   (2) the synthetic file    -> tune the simulator until the numbers and
%                                plots match the real ones
%
% WHAT IT PRODUCES
%   Console : a calibration summary (background rate, polarity balance,
%             target event count, tail shape, blur)
%   Fig. 1  : windowed snapshot  - net events (+1/-1) in a short time window
%   Fig. 2  : event rate vs time - noise floor and what the target adds
%   Fig. 3  : event offsets from the track - along-track and cross-track
%
% THE TRACK
% This script does not find the target; YOU tell it where the track is, as a
% straight line:   x(t) = trackX0 + trackVx*(t - trackT0)   (same for y).
%   * Synthetic file : copy the values the simulator prints at the end.
%   * Real file      : take them from your angular-rate pipeline's fit. The
%                      pipeline's px = [slope intercept] means
%                      trackVx = px(1), trackX0 = px(2), trackT0 = 0
%                      (same for py), IF your time zero matches this script's
%                      (see tZero_us below).

clear; clc; close all;

%% 1. INPUTS ---------------------------------------------------------------
fname     = 'synthetic_debris_events.txt';   % event file: x y polarity t_us
dataLabel = 'synthetic';                     % name shown in plot titles

nx = 1280;  ny = 720;       % sensor size (pixels)
zeroIndexed = false;        % true if the file's x,y start at 0 (Prophesee raw
                            % exports usually do); false if they start at 1

tZero_us = [];              % what counts as t = 0:
                            %   []  -> the file's first timestamp
                            %   0   -> keep the raw timestamps as they are
                            % Must match the time zero used for the track below.

% --- the track (straight line, pixels and seconds) ---
trackX0 = 1280;   trackY0 = 90;      % position at time trackT0 (px)
trackVx = -612.44; trackVy = 43.06;  % velocity (px/s)
trackT0 = 0;                         % time at which the track is at (X0,Y0) (s)

% --- plot / measurement settings ---
tc      = [];       % center time of the windowed snapshot (s). [] = middle of the file
winLen  = 5e-3;     % width of that window (s)
netClim = 10;       % snapshot color range (+/- this many net events per pixel)
zoomX   = 70;  zoomY = 35;   % snapshot zoom: half-width / half-height (px)
nTb     = 200;      % number of time bins for the rate plot

% "Gate" around the track used to count target events. Events inside it are
% treated as target-driven. Along-track range is negative BEHIND the target
% because the decay tail trails it and can be a few hundred pixels long.
alongRange = [-300 30];     % px along the track  (- = behind, + = ahead)
crossGate  = 15;            % px either side of the track
bgExclude  = 40;            % px: pixels closer than this to the track are
                            % ignored when measuring the background

%% 2. LOAD THE EVENTS ------------------------------------------------------
E   = readmatrix(fname);                 % columns: x y pol t_us
x   = E(:,1);  y = E(:,2);  pol = E(:,3);  tUs = E(:,4);
if zeroIndexed, x = x + 1;  y = y + 1; end      % convert to MATLAB 1-based
if min(pol) >= 0, pol = 2*pol - 1; end          % 0/1 polarity -> -1/+1

if max(x) > nx || max(y) > ny || min(x) < 1 || min(y) < 1
    error('Events fall outside a %dx%d sensor. Check nx, ny and zeroIndexed.', nx, ny);
end

if isempty(tZero_us), tZero_us = min(tUs); end
t = (tUs - tZero_us)*1e-6;                      % time in seconds

[t, ord] = sort(t);                             % make sure events are time-ordered
x = x(ord);  y = y(ord);  pol = pol(ord);
T   = max(t);                                   % duration of the recording (s)
nEv = numel(t);

if isempty(tc), tc = T/2; end

%% 3. GEOMETRY OF EVERY EVENT RELATIVE TO THE TRACK ------------------------
speed = hypot(trackVx, trackVy);
uHat  = [trackVx trackVy]/speed;                % unit vector along the motion

% Where the track says the target was at each event's time
xc = trackX0 + trackVx*(t - trackT0);
yc = trackY0 + trackVy*(t - trackT0);
dx = x - xc;   dy = y - yc;

% Split the offset into "along the motion" and "sideways"
along = dx*uHat(1) + dy*uHat(2);                % + = ahead of target, - = behind
cross = -dx*uHat(2) + dy*uHat(1);               % sideways distance from the path

% Events that count as target-driven
inGate = along >= alongRange(1) & along <= alongRange(2) & abs(cross) <= crossGate;
isPos  = pol > 0;
isNeg  = pol < 0;

% Perpendicular distance from the (infinite) track line, for the background
% measurement (independent of time)
distEv = abs((x - trackX0)*uHat(2) - (y - trackY0)*uHat(1));

%% 4. CALIBRATION SUMMARY --------------------------------------------------
fprintf('\n================ DIAGNOSTICS: %s ================\n', dataLabel);
fprintf('File: %s\n', fname);
fprintf('Events: %d   duration: %.3f s   overall rate: %.0f events/s\n', nEv, T, nEv/T);
fprintf('Overall polarity: %.1f%% positive\n', 100*mean(isPos));

% ---- Background (pixels away from the track) ----
Nimg = accumarray([y x], 1, [ny nx]);           % event count in every pixel
[X, Y] = meshgrid(1:nx, 1:ny);
distPix = abs((X - trackX0)*uHat(2) - (Y - trackY0)*uHat(1));
offMask = distPix > bgExclude;                  % pixels far from the track

rateOff   = Nimg(offMask)/T;                    % events/s for each background pixel
rs        = sort(rateOff);
meanBg    = mean(rateOff);
trimBg    = mean(rs(1:floor(0.99*numel(rs))));  % mean after dropping the busiest 1%
fracActive = mean(Nimg(offMask) > 0);           % fraction of pixels with >= 1 event
fracPosBg = mean(pol(distEv > bgExclude) > 0);  % + fraction of background events
fprintf('\n--- BACKGROUND (pixels > %d px from the track) ---\n', bgExclude);
fprintf('Mean rate            : %.4f events/s/pixel   <- use for noiseRate\n', meanBg);
fprintf('Mean, top 1%% removed : %.4f events/s/pixel   (ignores hot pixels / stars; reads LOW if noise is very sparse)\n', trimBg);
fprintf('Busiest pixel        : %.1f events/s\n', rs(end));
fprintf('Pixels with >=1 event: %.1f%%\n', 100*fracActive);
fprintf('Positive fraction    : %.1f%%\n', 100*fracPosBg);

% ---- Target events (inside the gate) ----
tt        = linspace(0, T, 10001)';             % when is the track on the sensor?
xt = trackX0 + trackVx*(tt - trackT0);   yt = trackY0 + trackVy*(tt - trackT0);
onSensor  = xt >= 1 & xt <= nx & yt >= 1 & yt <= ny;
transit   = sum(onSensor)*(T/10000);            % seconds the track is on the sensor
nGate     = nnz(inGate);
gateArea  = (alongRange(2) - alongRange(1))*(2*crossGate);
noiseInGate = meanBg*gateArea*T;                % noise events expected inside the gate
fprintf('\n--- TARGET (gate: along %+d..%+d px, cross +/-%d px) ---\n', ...
    alongRange(1), alongRange(2), crossGate);
fprintf('Track speed          : %.2f px/s\n', speed);
fprintf('Track on sensor for  : %.2f s\n', transit);
fprintf('Events in gate       : %d  (~%.0f events/s while on sensor)   <- compare with Phi_0 tuning\n', ...
    nGate, nGate/max(transit, eps));
fprintf('  of which noise est.: ~%.0f (%.1f%%)\n', noiseInGate, 100*noiseInGate/max(nGate,1));

for s = [1 -1]
    sel = inGate & (pol == s);
    if s > 0, name = '+'; else, name = '-'; end
    if ~any(sel), fprintf('%s events: none in gate\n', name); continue; end
    av = sort(along(sel));  nv = numel(av);
    q  = av(max(1, round([0.10 0.50 0.90]*nv)));            % 10th / 50th / 90th percentile
    fprintf('%s events: n = %6d | along-track median %+7.1f px (10%%/90%%: %+.1f / %+.1f) | rms cross-track %.2f px\n', ...
        name, nv, q(2), q(1), q(3), sqrt(mean(cross(sel).^2)));
end
fprintf('  (+ events sit at the target; cross-track rms grows with sigma_psf - compare real vs synthetic rather than reading it as sigma.)\n');
fprintf('  (- events trail behind; their along-track spread ~ tail length, set by a.)\n');

%% 5. FIGURE 1: WINDOWED SNAPSHOT ------------------------------------------
% Net event sum (+1/-1) per pixel inside a short time window. Red = net
% positive (target's leading edge), blue = net negative (decay tail).
cmap = [linspace(0,1,128)' linspace(0,1,128)' ones(128,1); ...
        ones(128,1) linspace(1,0,128)' linspace(1,0,128)'];   % blue-white-red

inWin   = t >= tc - winLen/2 & t < tc + winLen/2;
net_img = zeros(ny, nx);
if any(inWin)
    net_img = accumarray([y(inWin), x(inWin)], pol(inWin), [ny nx]);
else
    warning('No events in the snapshot window. Change tc.');
end
xcTc = trackX0 + trackVx*(tc - trackT0);        % where the track says the target is at tc
ycTc = trackY0 + trackVy*(tc - trackT0);

figure('Color','w','Name','Windowed snapshot','Position',[100 100 750 450]);
imagesc(net_img); axis image; colormap(gca, cmap);
try, clim([-netClim netClim]); catch, caxis([-netClim netClim]); end %#ok<CAXIS>
colorbar; hold on;
plot(xcTc, ycTc, 'k+', 'MarkerSize', 10);       % track position at time tc
xlim(xcTc + [-zoomX zoomX]);  ylim(ycTc + [-zoomY zoomY]);
xlabel('x [px]'); ylabel('y [px]');
title(sprintf('%s: net events in %.0f ms window at t = %.2f s', dataLabel, winLen*1e3, tc));

%% 6. FIGURE 2: EVENT RATE VS TIME -----------------------------------------
% Events per second, split by polarity. The flat floor is background noise;
% anything above it while the target is on the sensor is the target.
tb   = min(floor(t/T*nTb) + 1, nTb);            % time-bin index of each event
tAx  = ((1:nTb)-0.5)*T/nTb;                     % bin center times (s)
binW = T/nTb;
rPos = accumarray(tb(isPos), 1, [nTb 1])/binW;
rNeg = accumarray(tb(isNeg), 1, [nTb 1])/binW;

figure('Color','w','Name','Event rate vs time','Position',[100 100 800 400]);
plot(tAx, rPos, 'r-', tAx, rNeg, 'b-', tAx, rPos+rNeg, 'k-', 'LineWidth', 1);
grid on; xlabel('Time (s)'); ylabel('Event rate (events/s)');
legend('+ events', '- events', 'total', 'Location', 'best');
title(sprintf('%s: event rate vs time', dataLabel));

%% 7. FIGURE 3: EVENT OFFSETS FROM THE TRACK -------------------------------
% For events inside the gate: how far ahead (+) or behind (-) the track
% position each one landed. + events cluster at 0 (the target itself);
% - events trail behind it. The histogram width of the - events is the tail.
figure('Color','w','Name','Offsets from track','Position',[100 100 1150 450]);

subplot(1,2,1);
gi   = find(inGate);
sub  = gi(1:ceil(numel(gi)/20000):end);         % thin the points so plotting stays fast
posS = sub(isPos(sub));  negS = sub(isNeg(sub));
plot(t(negS), along(negS), 'b.', t(posS), along(posS), 'r.', 'MarkerSize', 4);
hold on; plot([0 T], [0 0], 'k--'); grid on;
xlabel('Time (s)'); ylabel('Along-track offset (px)  (+ = ahead)');
title(sprintf('%s: event offset from track', dataLabel));
legend('- events', '+ events', 'Location', 'best');

subplot(1,2,2);
histogram(along(inGate & isPos), 'BinWidth', 3, 'FaceColor', 'r', 'FaceAlpha', 0.5); hold on;
histogram(along(inGate & isNeg), 'BinWidth', 3, 'FaceColor', 'b', 'FaceAlpha', 0.5);
xlabel('Along-track offset (px)'); ylabel('Number of events'); grid on;
title('Where + and - events land'); legend('+ events', '- events');