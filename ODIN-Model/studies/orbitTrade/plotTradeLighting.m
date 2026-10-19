function summary = plotTradeLighting(res, nFirst)
% plotTradeLighting: sunlit, penumbra and umbra time of the trade cases, orbit by orbit and over the whole run.
% res: struct array from tradeLighting.m
% nFirst: number of orbits drawn in each timeline strip
% summary: struct array of mean full Sun, penumbra and umbra fractions, longest shadows and shadow dates per case

if nargin < 2
    nFirst = 5;
end

nCase = numel(res); % nCase: number of cases
col = [0.00 0.45 0.70
       0.84 0.37 0.00
       0.00 0.62 0.45
       0.80 0.47 0.65
       0.90 0.62 0.00
       0.34 0.71 0.91]; % col: 6 x 3 case colors matching the trade study figures
shade = [0.98 0.87 0.45
         0.60 0.60 0.66
         0.13 0.13 0.20]; % shade: 3 x 3 colors of full Sun, penumbra and umbra

fullSun = zeros(nCase, 1); % fullSun: mean fraction of time in full Sun [%]
penumbra = zeros(nCase, 1); % penumbra: mean fraction of time in penumbra only [%]
umbra = zeros(nCase, 1); % umbra: mean fraction of time in umbra [%]
longest = zeros(nCase, 1); % longest: longest shadow including penumbra [min]
longestUmbra = zeros(nCase, 1); % longestUmbra: longest umbra [min]
first = cell(nCase, 1); % first: date of the first shadow entry
last = cell(nCase, 1); % last: date of the last shadow exit

for k = 1:nCase
    [iv(k).p0, iv(k).p1] = pairIntervals(res(k).ecl.penumbraIn, res(k).ecl.penumbraOut, res(k).tEnd);
    [iv(k).u0, iv(k).u1] = pairIntervals(res(k).ecl.umbraIn, res(k).ecl.umbraOut, res(k).tEnd);
    % iv: 1 x nCase struct array of shadow (p0, p1) and umbra (u0, u1) interval start and end times [s]

    tShadow = sum(iv(k).p1 - iv(k).p0); % tShadow: total time in any shadow [s]
    tUmbra = sum(iv(k).u1 - iv(k).u0); % tUmbra: total time in umbra [s]

    fullSun(k) = 100 * (1 - tShadow/res(k).tEnd);
    umbra(k) = 100 * tUmbra/res(k).tEnd;
    penumbra(k) = 100 * (tShadow - tUmbra)/res(k).tEnd;
    longest(k) = max([0, iv(k).p1 - iv(k).p0]) / 60;
    longestUmbra(k) = max([0, iv(k).u1 - iv(k).u0]) / 60;

    if isempty(iv(k).p0)
        first{k} = '-';
        last{k} = '-';
    else
        first{k} = datestr(datenum(res(k).cfg.epoch) + iv(k).p0(1)/86400, 'dd mmm yyyy');
        last{k} = datestr(datenum(res(k).cfg.epoch) + iv(k).p1(end)/86400, 'dd mmm yyyy');
    end
end

summary = struct('name', {res.name}', 'fullSun', num2cell(fullSun), 'penumbra', num2cell(penumbra), ...
                 'umbra', num2cell(umbra), 'longest', num2cell(longest), 'longestUmbra', num2cell(longestUmbra), ...
                 'first', first, 'last', last);

fprintf('\ncase  full Sun [%%]  penumbra [%%]  umbra [%%]  longest shadow [min]  longest umbra [min]  first shadow  last shadow\n');
for k = 1:nCase
    fprintf('%-4s  %12.2f  %12.3f  %9.2f  %20.2f  %19.2f  %12s  %11s\n', res(k).name, fullSun(k), penumbra(k), ...
            umbra(k), longest(k), longestUmbra(k), first{k}, last{k});
end

figure('Color', 'w', 'Position', [80 80 1000 720]);

subplot(2, 1, 1);
tStart = zeros(1, nCase); % tStart: start of each strip, the run start [s]
drawStrips(res, iv, tStart, nFirst, shade, {res.name});
title(sprintf('First %d orbits from %s UTC, labels give full Sun per orbit', nFirst, ...
              datestr(datenum(res(1).cfg.epoch), 'dd mmm yyyy HH:MM')));
xlabel('orbits since epoch');

subplot(2, 1, 2);
lbl = cell(1, nCase); % lbl: row labels with the date of the longest shadow

for k = 1:nCase
    if isempty(iv(k).p0)
        tStart(k) = 0;
        lbl{k} = sprintf('%s  no shadow', res(k).name);
    else
        [~, m] = max(iv(k).p1 - iv(k).p0); % m: index of the longest shadow
        tMid = (iv(k).p0(m) + iv(k).p1(m)) / 2; % tMid: centre of the longest shadow [s]

        % Strip start placing the longest shadow at the centre of the middle orbit, kept inside the run:
        tStart(k) = min(max(0, tMid - (floor(nFirst/2) + 0.5)*res(k).T), max(0, res(k).tEnd - nFirst*res(k).T));

        lbl{k} = sprintf('%s  %s', res(k).name, datestr(datenum(res(k).cfg.epoch) + tMid/86400, 'dd mmm yy'));
    end
end

drawStrips(res, iv, tStart, nFirst, shade, lbl);
title('Orbits around the longest shadow of each case');
xlabel('orbits from the start of the strip');

figure('Color', 'w', 'Position', [120 60 1000 720]);

for k = 1:nCase
    subplot(ceil(nCase/2), 2, k);

    d = datenum(res(k).cfg.epoch) + iv(k).p0/86400; % d: 1 x M shadow entry dates [days]
    dur = (iv(k).p1 - iv(k).p0) / 60; % dur: 1 x M shadow durations including penumbra [min]

    % Season break where consecutive shadows start more than two orbits apart:
    gap = [false, diff(iv(k).p0) > 2*res(k).T]; % gap: 1 x M mask of shadows that start a new season
    d(gap) = NaN;
    dur(gap) = NaN;

    c = double(res(k).name) - double('A') + 1; % c: row of col for this case letter
    plot(d, dur, 'Color', col(c,:), 'LineWidth', 1.5); grid on;
    xlim(datenum(res(k).cfg.epoch) + [0 res(k).tEnd/86400]);
    ylim([0 40]);
    datetick('x', 'mmm yy', 'keeplimits');
    ylabel('shadow [min]');
    title(sprintf('%s   full Sun %.1f%%   longest %.1f min', res(k).label, fullSun(k), longest(k)), ...
          'FontWeight', 'normal', 'FontSize', 9);
end

end


function [t0, t1] = pairIntervals(tIn, tOut, tEnd)
% pairIntervals: entry and exit times paired into intervals, closed at the run start and end.
% tIn: 1 x K entry times [s]
% tOut: 1 x L exit times [s]
% tEnd: run length [s]
% t0: 1 x M interval start times [s]
% t1: 1 x M interval end times [s]

t0 = tIn;
t1 = tOut;

if ~isempty(t1) && (isempty(t0) || t1(1) < t0(1))
    t0 = [0, t0];
end

if numel(t0) > numel(t1)
    t1 = [t1, tEnd];
end

end


function drawStrips(res, iv, tStart, nOrb, shade, lbl)
% drawStrips: one timeline row per case coloured by full Sun, penumbra and umbra, labelled with full Sun per orbit.
% res: struct array from tradeLighting.m
% iv: struct array of shadow and umbra intervals [s]
% tStart: 1 x nCase strip start times [s]
% nOrb: number of orbits in each strip
% shade: 3 x 3 colors of full Sun, penumbra and umbra
% lbl: 1 x nCase row labels

nCase = numel(res); % nCase: number of cases
hold on;

for k = 1:nCase
    y = nCase - k + 1; % y: row centre
    T = res(k).T; % T: orbital period [s]
    tEnd = min(tStart(k) + nOrb*T, res(k).tEnd); % tEnd: strip end time, clipped to the run [s]

    rect(0, (tEnd - tStart(k))/T, y, shade(1,:));

    for q = find(iv(k).p1 > tStart(k) & iv(k).p0 < tEnd)
        rect((max(iv(k).p0(q), tStart(k)) - tStart(k))/T, (min(iv(k).p1(q), tEnd) - tStart(k))/T, y, shade(2,:));
    end

    for q = find(iv(k).u1 > tStart(k) & iv(k).u0 < tEnd)
        rect((max(iv(k).u0(q), tStart(k)) - tStart(k))/T, (min(iv(k).u1(q), tEnd) - tStart(k))/T, y, shade(3,:));
    end

    for j = 1:floor((tEnd - tStart(k))/T + 1e-9)
        w0 = tStart(k) + (j-1)*T; % w0: orbit window start [s]
        w1 = w0 + T; % w1: orbit window end [s]

        % Interval overlap with the orbit window, max(0, min(t1, w1) - max(t0, w0)):
        s = sum(max(0, min(iv(k).p1, w1) - max(iv(k).p0, w0))); % s: shadow time in this orbit [s]

        text(j - 0.5, y + 0.3, sprintf('%.1f%%', 100*(1 - s/T)), ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 8);
    end
end

for j = 1:nOrb-1
    plot([j j], [0.5 nCase+0.5], ':', 'Color', [0.4 0.4 0.4]);
end

h = []; % h: legend handles of full Sun, penumbra and umbra
for c = 1:3
    h(c) = patch(NaN, NaN, shade(c,:), 'EdgeColor', 'none');
end
legend(h, {'full Sun', 'penumbra', 'umbra'}, 'Location', 'eastoutside');

xlim([0 nOrb]);
ylim([0.5 nCase + 0.65]);
set(gca, 'YTick', 1:nCase, 'YTickLabel', fliplr(lbl), 'Box', 'on');

end


function rect(x0, x1, y, c)
% rect: filled horizontal bar of one strip row.
% x0: bar start [orbits]
% x1: bar end [orbits]
% y: row centre
% c: 1 x 3 fill color

patch([x0 x1 x1 x0], y + 0.25*[-1 -1 1 1], c, 'EdgeColor', 'none');

end
