% eclipseCrossCheck: toolbox eclipse intervals against eclipseTimes for trade case B at 60 s and 1 s scenario sample time.

here = fileparts(mfilename('fullpath')); % here: folder holding this script
addpath(fullfile(here, '..'));
odinPath;

outDir = fullfile(here, 'results'); % outDir: folder for logs and figures
if ~isfolder(outDir)
    mkdir(outDir);
end
logFile = fullfile(outDir, 'eclipseCrossCheck.txt'); % logFile: text log of this run
if isfile(logFile)
    delete(logFile);
end
diary(logFile);

nDay = 2; % nDay: run length [days]
res = tradeLighting(nDay, 2); % res: case B parameters and shadow times refined from the 60 s grid
cfg = res.cfg; % cfg: parameter struct of case B
cfg.ode.dtOut = 1;

out = propagate(cfg, nDay*86400/res.T); % out: state history of case B on a 1 s grid
ecl = eclipseTimes(out, cfg); % ecl: shadow times refined from the 1 s grid [s]

[p0, p1] = pairShadow(ecl.penumbraIn, ecl.penumbraOut); % p0, p1: 1 x K complete shadow interval starts and ends [s]
[q0, q1] = pairShadow(res.ecl.penumbraIn, res.ecl.penumbraOut); % q0, q1: the same from the 60 s grid [s]
fprintf('\n%s, %d complete shadows\n', res.label, numel(p0));
fprintf('eclipseTimes from the 60 s and 1 s grids, largest entry and exit difference: %.2e s, %.2e s\n', ...
        max(abs(q0 - p0)), max(abs(q1 - p1)));
fprintf('shadow duration %.2f to %.2f min\n', min(p1 - p0)/60, max(p1 - p0)/60);

tt = datetime(cfg.epoch) + seconds(out.t(:)); % tt: N x 1 UTC sample times
dts = [60 1]; % dts: scenario sample times [s]
nm = {'corrected', 'uncorrected'}; % nm: satellite names, rotation applied and omitted
dS = cell(2, 2); % dS: toolbox minus eclipseTimes shadow start per eclipse, by sample time and satellite [s]
dE = cell(2, 2); % dE: toolbox minus eclipseTimes shadow end per eclipse, by sample time and satellite [s]
dD = cell(2, 2); % dD: toolbox minus eclipseTimes shadow duration per eclipse, by sample time and satellite [s]

for j = 1:2
    sc = satelliteScenario(tt(1), tt(end), dts(j)); % sc: scenario at this sample time
    sat = [odinToScenario(sc, out, nm{1}), ...
           satellite(sc, timetable(tt, 1e3*out.r.'), timetable(tt, 1e3*out.v.'), 'CoordinateFrame', 'inertial', 'Name', nm{2})];
    % sat: 1 x 2 satellites with the rotation applied and omitted

    for s = 1:2
        e = eclipse(sat(s)); % e: eclipse analysis object with the default dual-cone model
        iv = eclipseIntervals(e); % iv: toolbox shadow intervals

        s0 = seconds(iv.StartTime - sc.StartTime).'; % s0: 1 x L toolbox shadow starts since epoch [s]
        s1 = seconds(iv.EndTime - sc.StartTime).'; % s1: 1 x L toolbox shadow ends since epoch [s]
        D = iv.Duration.'; % D: 1 x L toolbox shadow durations [s]
        if isduration(D)
            D = seconds(D);
        end

        ok = p0 > dts(j) & p1 < out.t(end) - dts(j); % ok: 1 x K mask of shadows clear of the scenario ends
        m = matchIntervals(p0(ok), p1(ok), s0, s1); % m: index of the toolbox interval overlapping each shadow

        dS{j,s} = s0(m) - p0(ok);
        dE{j,s} = s1(m) - p1(ok);
        dD{j,s} = D(m) - (p1(ok) - p0(ok));

        fprintf('\n%2d s sample time, %s: %d of %d shadows matched\n', dts(j), nm{s}, numel(m), nnz(ok));
        fprintf('  start  mean %+7.2f s  min %+7.2f  max %+7.2f\n', mean(dS{j,s}), min(dS{j,s}), max(dS{j,s}));
        fprintf('  end    mean %+7.2f s  min %+7.2f  max %+7.2f\n', mean(dE{j,s}), min(dE{j,s}), max(dE{j,s}));
        fprintf('  dur    mean %+7.2f s  min %+7.2f  max %+7.2f\n', mean(dD{j,s}), min(dD{j,s}), max(dD{j,s}));
        fprintf('  largest offset of a toolbox start from the sample grid %.2e s, of a duration from a whole number of samples %.2e s\n', ...
                max(abs(s0 - dts(j)*round(s0/dts(j)))), max(abs(D - dts(j)*round(D/dts(j)))));

        if dts(j) == 1
            f = eclipseStatus(e); % f: 1 x N toolbox visible fraction of the solar disk
            g = gradient(out.illum, out.t); % g: 1 x N rate of change of the eclipseTimes model fraction [1/s]
            w = out.illum > 0.02 & out.illum < 0.98; % w: 1 x N mask of penumbra samples

            % First order time shift of the toolbox fraction against the model fraction:
            lag = (out.illum(w) - f(w)) ./ g(w); % lag: toolbox transition lag behind the model at each penumbra sample [s]

            fprintf('  toolbox lag from the penumbra fraction: entry mean %+6.3f s (sd %.3f), exit mean %+6.3f s (sd %.3f)\n', ...
                    mean(lag(g(w) < 0)), std(lag(g(w) < 0)), mean(lag(g(w) > 0)), std(lag(g(w) > 0)));
            fprintf('  largest fraction difference %.4f\n', max(abs(out.illum - f)));
        end
    end
end

figure('Color', 'w', 'Position', [100 100 900 600]);

subplot(2, 1, 1);
plot(dD{1,1}, 'o-', 'LineWidth', 1.2); grid on; hold on;
plot(dD{2,1}, 's-', 'LineWidth', 1.2);
ylabel('toolbox - eclipseTimes [s]'); xlabel('shadow number');
legend('duration, 60 s sample time', 'duration, 1 s sample time', 'Location', 'best');
title(sprintf('%s, rotation applied', res.label), 'FontWeight', 'normal');

subplot(2, 1, 2);
plot(dS{2,1}, 'o-', 'LineWidth', 1.2); grid on; hold on;
plot(dE{2,1}, 's-', 'LineWidth', 1.2);
plot(dS{2,2}, 'o--', 'LineWidth', 1.2);
plot(dE{2,2}, 's--', 'LineWidth', 1.2);
ylabel('toolbox - eclipseTimes [s]'); xlabel('shadow number');
legend('start, rotation applied', 'end, rotation applied', 'start, rotation omitted', 'end, rotation omitted', 'Location', 'best');
title('1 s sample time', 'FontWeight', 'normal');

exportgraphics(gcf, fullfile(outDir, 'eclipseCrossCheck.png'), 'Resolution', 150);

diary off;
save(fullfile(outDir, 'eclipseCrossCheck.mat'), 'dS', 'dE', 'dD', 'dts', 'nm', 'p0', 'p1');


function [t0, t1] = pairShadow(tIn, tOut)
% pairShadow: each entry paired with the exit that follows it, incomplete intervals dropped.
% tIn: 1 x K entry times [s]
% tOut: 1 x L exit times [s]
% t0: 1 x M interval start times [s]
% t1: 1 x M interval end times [s]

t0 = zeros(1, 0);
t1 = zeros(1, 0);

for j = 1:numel(tIn)
    n = find(tOut > tIn(j), 1); % n: index of the first exit after this entry
    if ~isempty(n)
        t0(end+1) = tIn(j);
        t1(end+1) = tOut(n);
    end
end

end


function m = matchIntervals(p0, p1, s0, s1)
% matchIntervals: toolbox interval with the largest overlap for each reference interval.
% p0: 1 x K reference interval starts [s]
% p1: 1 x K reference interval ends [s]
% s0: 1 x L toolbox interval starts [s]
% s1: 1 x L toolbox interval ends [s]
% m: 1 x K index of the matching toolbox interval

m = zeros(1, numel(p0));

for j = 1:numel(p0)
    % Interval overlap, min(ends) - max(starts):
    ov = min(s1, p1(j)) - max(s0, p0(j)); % ov: 1 x L overlap with each toolbox interval [s]
    [~, m(j)] = max(ov);
end

end
