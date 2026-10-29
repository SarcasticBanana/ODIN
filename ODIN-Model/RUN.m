% RUN: the ODIN constellation over nDay days, orbit, shadow, array power, ground track and downlink figures, the contact table, downlink times in seconds, then the 3D viewer.

clear; close all; clc;

here = fileparts(mfilename('fullpath')); % here: repository root folder
addpath(here);
odinPath;

nDay = 7; % nDay: run length [days]
nView = 24; % nView: 3D viewer scenario length [h]
nOrbP = 10; % nOrbP: orbits of ODIN 1.1 given to arrayPower
capMin = 120; % capMin: downlink time per satellite per day in the battery trade [min]
outDir = fullfile(here, 'results'); % outDir: folder for the contact table
if ~isfolder(outDir)
    mkdir(outDir);
end

col = [0.00 0.45 0.70
       0.84 0.37 0.00
       0.00 0.62 0.45
       0.80 0.47 0.65
       0.90 0.62 0.00
       0.34 0.71 0.91]; % col: 6 x 3 plane colours, as in scenarioDemo.m
gsCol = col([5 4 6], :); % gsCol: 3 x 3 ground station colours

base = config(); % base: parameter struct of ODIN 1.1 from the Inputs sheet
st = stations(); % st: ground stations from the Stations sheet
[cfgs, names] = constellation(base);
% cfgs: parameter struct of each satellite from constellation.m
% names: name of each satellite

nSat = numel(cfgs); % nSat: number of satellites
nG = numel(st); % nG: number of ground stations
nP = base.const.nPerPlane; % nP: satellites per plane
w = 0.6*(0:nP-1)/max(nP - 1, 1); % w: share of white in the shade of each slot
t0 = datetime(base.epoch); % t0: UTC instant at sim start
tag = sprintf('%g days from %s UTC', nDay, datestr(datenum(base.epoch), 'yyyy-mm-dd HH:MM')); % tag: run span for the figure titles

out = cell(1, nSat); % out: state history of each satellite
ecl = cell(1, nSat); % ecl: shadow entry and exit times of each satellite
con = cell(nSat, nG); % con: contact windows of each satellite with each station
dl = cell(1, nSat); % dl: downlink start and end times of each satellite from downlinkTimes.m, K x 2 [s]
c = zeros(nSat, 3); % c: colour of each satellite
for j = 1:nSat
    k = ceil(j/nP); % k: plane of satellite j
    s = j - (k - 1)*nP; % s: slot of satellite j in its plane

    % Plane colour mixed with white for the slot shade:
    c(j,:) = (1 - w(s))*col(k,:) + w(s);

    % Kepler's third law:
    T = 2*pi*sqrt((cfgs(j).Re + cfgs(j).h)^3/cfgs(j).mu); % T: orbital period of satellite j [s]

    out{j} = propagate(cfgs(j), (nDay*86400 + 1)/T);
    if base.use.eclipse
        ecl{j} = eclipseTimes(out{j}, cfgs(j));
    end
    [dl{j}, con(j, :)] = downlinkTimes(out{j}, cfgs(j), st);
    fprintf('%s propagated, %d contact windows\n', names{j}, sum(cellfun(@(x) numel(x.aos), con(j, :))));
end

tEnd = out{1}.t(end); % tEnd: last output time [s]
ps = plotScale(out{1}, base); % ps: axis units and averaging settings for this run length
jd0 = datenum(base.epoch) + 1721058.5; % jd0: Julian date at sim start [days]

figure('Color', 'w', 'Name', 'Orbits', 'Position', [60 60 1000 650]);
ax = subplot(2, 1, 1); % ax: axes being drawn
hold(ax, 'on'); grid(ax, 'on');
hO = gobjects(1, nSat); % hO: altitude line of each satellite, for the legend
for j = 1:nSat
    alt = vecnorm(out{j}.r) - base.Re; % alt: 1 x N altitude above cfg.Re [km]
    if ps.dense
        altLo = movmin(alt, ps.wOrb); % altLo: 1 x N per orbit minimum altitude [km]
        altHi = movmax(alt, ps.wOrb); % altHi: 1 x N per orbit maximum altitude [km]
        hO(j) = plot(ax, ps.t(ps.dec), altHi(ps.dec), 'Color', c(j,:), 'LineWidth', 1.0);
        plot(ax, ps.t(ps.dec), altLo(ps.dec), 'Color', c(j,:), 'LineWidth', 1.0);
    else
        hO(j) = plot(ax, ps.t, alt, 'Color', c(j,:), 'LineWidth', 1.0);
    end
end
ylabel(ax, 'altitude above cfg.Re [km]');
if ps.dense
    title(ax, sprintf('per orbit minimum and maximum altitude, %s', tag));
else
    title(ax, sprintf('altitude, %s', tag));
end
legend(ax, hO, names, 'Location', 'eastoutside');

ax = subplot(2, 1, 2);
hold(ax, 'on'); grid(ax, 'on');
for j = 1:nSat
    raan = rad2deg(atan2(out{j}.hvec(1,:), -out{j}.hvec(2,:))); % raan: 1 x N right ascension of the ascending node [deg]

    % Mean Sun right ascension as in tradeLighting.m, L = 280.460 + 36000.771 T deg:
    raSun = 280.460 + 36000.771*(jd0 + out{j}.t/86400 - 2451545.0)/36525; % raSun: 1 x N mean Sun right ascension [deg]

    % Local time of the ascending node, LTAN = 12 h + (RAAN - RA_sun)/(15 deg/h):
    ltan = mod(12 + (raan - raSun)/15, 24); % ltan: 1 x N local time of the ascending node [h]

    plot(ax, ps.t(ps.dec), ltan(ps.dec), 'Color', c(j,:), 'LineWidth', 1.0);
end
ylabel(ax, 'LTAN [h]'); xlabel(ax, sprintf('time [%s]', ps.tName));
title(ax, 'local time of the ascending node, mean Sun');
legend(ax, names, 'Location', 'eastoutside');

if base.use.eclipse
    figure('Color', 'w', 'Name', 'Shadow', 'Position', [80 60 1000 750]);
    ax = subplot(3, 1, 1);
    hold(ax, 'on'); grid(ax, 'on');
    worst = [0 0 0]; % worst: satellite, entry time [s] and length [s] of the longest shadow
    for j = 1:nSat
        d = ecl{j}.penumbraDur; % d: 1 x K shadow lengths from penumbra entry to exit [s]
        plot(ax, ecl{j}.penumbraIn/ps.tDiv, d/60, 'o-', 'Color', c(j,:), 'MarkerFaceColor', c(j,:), 'MarkerSize', 4);
        [dMax, iMax] = max([0, d]);
        % dMax: longest shadow of satellite j, zero without shadow [s]
        % iMax: index of that shadow counted after the leading zero
        if dMax > worst(3)
            worst = [j, ecl{j}.penumbraIn(iMax - 1), dMax];
        end
    end
    if worst(1) == 0
        text(ax, 0.5, 0.5, 'no Earth shadow in this run', 'Units', 'normalized', 'HorizontalAlignment', 'center');
    end
    xlim(ax, [0 ps.t(end)]);
    ylabel(ax, 'shadow per orbit [min]');
    title(ax, sprintf('Earth shadow, penumbra entry to exit, %s', tag));

    ax = subplot(3, 1, 2);
    hold(ax, 'on'); grid(ax, 'on');
    for j = 1:nSat
        % Sun elevation above the orbit plane, beta = asin(h_hat . s_hat):
        beta = asind(sum(out{j}.hvec./vecnorm(out{j}.hvec).*out{j}.sunhat, 1)); % beta: 1 x N beta angle [deg]

        plot(ax, ps.t(ps.dec), beta(ps.dec), 'Color', c(j,:), 'LineWidth', 1.0);
    end

    % Shadow-free beta angle of a circular orbit, |beta| > asin(Re/a):
    bStar = asind(base.Re/(base.Re + base.h)); % bStar: beta angle beyond which an orbit at cfg.h sees no shadow [deg]

    plot(ax, ps.t([1 end]), [bStar bStar], 'k--', ps.t([1 end]), -[bStar bStar], 'k--');
    ylabel(ax, 'beta angle [deg]');
    title(ax, 'Sun angle above the orbit plane, dashed: no shadow outside this band');

    ax = subplot(3, 1, 3);
    hold(ax, 'on'); grid(ax, 'on');
    for j = 1:nSat
        S = movmean(solarFlux(out{j}, cfgs(j)), ps.wOrb); % S: 1 x N one orbit running mean irradiance at the spacecraft [W/m^2]
        plot(ax, ps.t(ps.dec), S(ps.dec), 'Color', c(j,:), 'LineWidth', 1.0);
    end
    ylabel(ax, 'orbit average [W/m^2]'); xlabel(ax, sprintf('time [%s]', ps.tName));
    title(ax, 'irradiance on a surface held normal to the Sun');

    if worst(1) == 0
        fprintf('\nno Earth shadow for any satellite over the %g days\n', nDay);
    else
        fprintf('\nlongest shadow: %s, %.2f min from %s UTC\n', names{worst(1)}, worst(3)/60, ...
                char(string(t0 + seconds(worst(2)), 'yyyy-MM-dd HH:mm:ss')));
    end

    outP = propagate(base, nOrbP); % outP: state history of ODIN 1.1 over nOrbP orbits
    ar = arrayPower(outP, base, base.power.panelArea, struct('incidenceSweepArea_m2', base.power.panelArea, 'makePlots', true));
    % ar: area study and incidence sweep at cfg.power.panelArea from arrayPower.m
    fprintf('array at %g m^2, %d orbits of %s, normal Sun: peak %.2f W, window average %.2f W at the bus, ageing allowance included\n', ...
            base.power.panelArea, nOrbP, names{1}, ar.summary.PeakPlanning_W, ar.summary.WindowAvgPlanning_W);
end

figure('Color', 'w', 'Name', 'Ground track', 'Position', [100 60 1100 620]);
gx = geoaxes; % gx: map axes
geobasemap(gx, 'darkwater');
hold(gx, 'on');
e2 = base.f*(2 - base.f); % e2: squared eccentricity of the ellipsoid
bE = base.Re*(1 - base.f); % bE: polar radius of the ellipsoid [km]
hM = gobjects(1, nSat + nG); % hM: line of each satellite and station, for the legend
for j = 1:nSat
    q1 = out{j}.t <= out{j}.T; % q1: samples of the first orbit
    r = out{j}.r(:, q1); % r: 3 x M position over the first orbit [km]
    th = earthRotation(out{j}.t(q1), base); % th: 1 x M Greenwich angle [rad]

    % Model frame to Earth fixed, rotation about z by the Greenwich angle:
    xe = cos(th).*r(1,:) + sin(th).*r(2,:); % xe: 1 x M Earth fixed x [km]
    ye = cos(th).*r(2,:) - sin(th).*r(1,:); % ye: 1 x M Earth fixed y [km]

    p = hypot(xe, ye); % p: 1 x M distance from the polar axis [km]
    u = atan2(base.Re*r(3,:), bE*p); % u: 1 x M parametric latitude [rad]

    % Geodetic latitude, Bowring's formula:
    lat = atan2d(r(3,:) + e2/(1 - e2)*bE*sin(u).^3, p - e2*base.Re*cos(u).^3); % lat: 1 x M geodetic latitude of the subsatellite point [deg]

    lon = atan2d(ye, xe); % lon: 1 x M longitude of the subsatellite point, east positive [deg]

    kw = find(abs(diff(lon)) > 180); % kw: samples followed by a wrap of the longitude
    [~, ord] = sort([1:numel(lon), kw + 0.5]); % ord: sample order with a break after each wrap
    lat = [lat, NaN(size(kw))];
    lon = [lon, NaN(size(kw))];
    hM(j) = geoplot(gx, lat(ord), lon(ord), '-', 'Color', c(j,:), 'LineWidth', 1.2);
    geoplot(gx, lat(1), lon(1), 'o', 'Color', c(j,:), 'MarkerFaceColor', c(j,:), 'MarkerSize', 7);
end

az = linspace(0, 360, 361); % az: azimuths around each station [deg]
for g = 1:nG
    cg = gsCol(mod(g - 1, size(gsCol, 1)) + 1, :); % cg: colour of station g

    % Earth central angle from the station to the edge of view at its mask, satellite at cfg.h:
    psi = 90 - st(g).mask - asind(base.Re*cosd(st(g).mask)/(base.Re + base.h)); % psi: radius of the visibility circle [deg]

    % Small circle of radius psi about the station:
    latC = asind(sind(st(g).lat)*cosd(psi) + cosd(st(g).lat)*sind(psi)*cosd(az)); % latC: 1 x 361 circle latitudes [deg]
    lonC = mod(st(g).lon + atan2d(sind(az)*sind(psi)*cosd(st(g).lat), cosd(psi) - sind(st(g).lat)*sind(latC)) + 180, 360) - 180;
    % lonC: 1 x 361 circle longitudes [deg]

    kw = find(abs(diff(lonC)) > 180);
    [~, ord] = sort([1:numel(lonC), kw + 0.5]);
    latC = [latC, NaN(size(kw))];
    lonC = [lonC, NaN(size(kw))];
    hM(nSat + g) = geoplot(gx, latC(ord), lonC(ord), '--', 'Color', cg, 'LineWidth', 1.2);
    geoplot(gx, st(g).lat, st(g).lon, '^', 'Color', cg, 'MarkerFaceColor', cg, 'MarkerSize', 9);
end
legend(gx, hM, [names, {st.name}], 'Location', 'eastoutside');
title(gx, sprintf('first orbit of each satellite from %s UTC, dashed: edge of view at each station''s mask, %g km', ...
                  datestr(datenum(base.epoch), 'yyyy-mm-dd HH:MM'), base.h));

figure('Color', 'w', 'Name', 'Downlink', 'Position', [120 60 1100 700]);
ax = subplot(3, 1, [1 2]);
hold(ax, 'on'); grid(ax, 'on');
hL = gobjects(1, nG); % hL: line of each station, for the legend
for g = 1:nG
    cg = gsCol(mod(g - 1, size(gsCol, 1)) + 1, :);
    for j = 1:nSat
        a = max(con{j, g}.aos, 0); % a: 1 x K window starts, the run start for a window open at the start [s]
        b = min(con{j, g}.los, tEnd); % b: 1 x K window ends, the run end for a window open at the end [s]
        hL(g) = plot(ax, t0 + seconds(reshape([a; b; NaN(size(a))], 1, [])), reshape([j; j; NaN]*ones(size(a)), 1, []), ...
                     '-', 'Color', cg, 'LineWidth', 7);
    end
end
set(ax, 'YDir', 'reverse', 'YTick', 1:nSat, 'YTickLabel', names, 'YLim', [0.4, nSat + 0.6]);
legend(ax, hL, {st.name}, 'Location', 'northoutside', 'Orientation', 'horizontal');
title(ax, sprintf('windows above each station''s mask, %s, times to the millisecond in results/contacts.csv', tag));

ax = subplot(3, 1, 3);
hold(ax, 'on'); grid(ax, 'on');
nD = ceil(tEnd/86400); % nD: days in the run, the last one partial when the run is not whole days
d0 = (0:nD-1)*86400; % d0: 1 x nD day start times [s]
dayMin = zeros(nD, nSat); % dayMin: nD x nSat time in view of any station per day [min]
for j = 1:nSat
    % Overlap of each downlink window with each day, max(0, min(ends) - max(starts)):
    dayMin(:, j) = sum(max(0, min(dl{j}(:, 2), d0 + 86400) - max(dl{j}(:, 1), d0)), 1).'/60;
end
hb = bar(ax, 1:nD, dayMin, 'grouped'); % hb: bar series of each satellite
for j = 1:nSat
    set(hb(j), 'FaceColor', c(j,:), 'EdgeColor', 'none');
end
plot(ax, [0.5, nD + 0.5], [capMin capMin], 'k--', 'LineWidth', 1.2);
set(ax, 'XTick', 1:nD, 'XTickLabel', cellstr(datestr(datenum(base.epoch) + (0:nD-1), 'dd mmm')), 'XLim', [0.5, nD + 0.5]);
ylabel(ax, 'in view per day [min]');
title(ax, sprintf('time in view of any station per day from sim start, dashed: %g min per day of downlink in the battery trade', capMin));

W = zeros(0, 6); % W: one row per window, satellite, station, acquisition [s], loss [s], highest elevation [deg] and its time [s]
for j = 1:nSat
    for g = 1:nG
        K = numel(con{j, g}.aos); % K: number of windows of satellite j with station g
        W = [W; repmat([j g], K, 1), con{j, g}.aos(:), con{j, g}.los(:), con{j, g}.elMax(:), con{j, g}.tMax(:)];
    end
end
key = W(:, 3); % key: acquisition times with windows open at the start first [s]
key(isnan(key)) = -Inf;
[~, iw] = sort(key); % iw: windows in order of acquisition
W = W(iw, :);

utc = @(x) string(t0 + milliseconds(round(1000*x)), 'yyyy-MM-dd HH:mm:ss.SSS'); % utc: UTC text of times since sim start, to the millisecond
C = table(names(W(:, 1)).', {st(W(:, 2)).name}.', utc(W(:, 3)), utc(W(:, 4)), W(:, 4) - W(:, 3), W(:, 5), utc(W(:, 6)), ...
          'VariableNames', {'satellite', 'station', 'AOS_UTC', 'LOS_UTC', 'duration_s', 'maxElevation_deg', 'maxElevation_UTC'});
% C: contact table, one row per window, empty times where a window runs past either end of the run
csvFile = fullfile(outDir, 'contacts.csv'); % csvFile: contact table file
writetable(C, csvFile);

fprintf('\n%d contact windows over %s, all of them in %s\n', size(W, 1), tag, csvFile);
fprintf('downlink start and end times of satellite j in seconds from sim start: dl{j}, one row per window, any station\n');
fprintf('first complete windows:\n');
for q = find(~isnan(W(:, 3)) & ~isnan(W(:, 4)), 5).'
    fprintf('  %-9s %-16s AOS %s  LOS %s  %8.3f s  max elevation %5.2f deg\n', names{W(q, 1)}, st(W(q, 2)).name, ...
            char(utc(W(q, 3))), char(utc(W(q, 4))), W(q, 4) - W(q, 3), W(q, 5));
end

fprintf('\nper day over the run\n%-10s', '');
for g = 1:nG
    fprintf('  %-18s', st(g).name);
end
fprintf('  %12s  days over\n%-10s', 'any station', 'satellite');
for g = 1:nG
    fprintf('  %8s  %8s', 'windows', 'minutes');
end
fprintf('  %12s  %g min\n', 'minutes', capMin);
for j = 1:nSat
    fprintf('%-10s', names{j});
    for g = 1:nG
        dur = min(con{j, g}.los, tEnd) - max(con{j, g}.aos, 0); % dur: 1 x K window lengths inside the run [s]
        fprintf('  %8.2f  %8.1f', numel(dur)/(tEnd/86400), sum(dur)/60/(tEnd/86400));
    end
    fprintf('  %12.1f  %d of %d\n', sum(dayMin(:, j))/(tEnd/86400), nnz(dayMin(:, j) > capMin), nD);
end

scenarioDemo([], [], nView);
