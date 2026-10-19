% numericalCrossCheck: J2 model against the toolbox numerical propagator for the six trade cases, toolbox runs in parallel.

here = fileparts(which('numericalCrossCheck')); % here: folder holding this script
addpath(fullfile(here, '..'));
odinPath;

outDir = fullfile(here, 'results'); % outDir: folder for logs and figures
if ~isfolder(outDir)
    mkdir(outDir);
end

if ~exist('nDay', 'var')
    nDay = 365; % nDay: run length [days]
end

logFile = fullfile(outDir, sprintf('numericalCrossCheck_%dd.txt', nDay)); % logFile: text log of this run
if isfile(logFile)
    delete(logFile);
end
diary(logFile);

nCase = 6; % nCase: number of trade cases
deg = 120; % deg: EGM2008 degree and order, the toolbox default
tol = 1e-8; % tol: relative and absolute tolerance of the toolbox ODE solver, the toolbox default
swFile = fullfile(outDir, 'SW-All.csv'); % swFile: CelesTrak consolidated space weather file, observed and predicted

dt = 60; % dt: output spacing [s]
t = 0 : dt : nDay*86400; % t: 1 x N output times since epoch [s]
N = numel(t); % N: number of output samples

ckDir = fullfile(tempdir, 'odin_numericalCrossCheck'); % ckDir: checkpoint folder, so an interrupted run resumes where it stopped
if ~isfolder(ckDir)
    mkdir(ckDir);
end

refFile = fullfile(ckDir, sprintf('J2model_%dd.mat', nDay)); % refFile: checkpoint of the J2 model runs
if isfile(refFile)
    load(refFile, 'ref', 'R', 'x0', 'altJ2', 'lowJ2');
    fprintf('J2 model runs loaded from %s\n', refFile);
else
    ref = struct('name', {}, 'label', {}, 'cfg', {}, 'T', {}, 'tEnd', {}, 'ecl', {}); % ref: J2 model runs in the tradeLighting.m result format
    R = zeros(3, 3, nCase); % R: mean equator and equinox of epoch to GCRF rotation of each case
    x0 = zeros(6, nCase); % x0: initial GCRF state of each case [m; m/s]
    altJ2 = zeros(N, nCase); % altJ2: N x nCase orbit averaged altitude of the J2 model [km]
    lowJ2 = zeros(1, nCase); % lowJ2: lowest altitude of the J2 model over the run [km]

    for k = 1:nCase
        r1 = tradeLighting(1/24, k); % r1: case parameters from tradeLighting.m, one hour run
        cfg = r1.cfg; % cfg: parameter struct of this case

        tic;
        out = propagate(cfg, (nDay*86400 + 1)/r1.T); % out: J2 model state history, one second past the last output time
        out.t = t;
        out.r = out.r(:, 1:N);
        out.rsun = out.rsun(:, 1:N);

        ref(k) = struct('name', r1.name, 'label', r1.label, 'cfg', cfg, 'T', out.T, 'tEnd', t(end), 'ecl', eclipseTimes(out, cfg));
        fprintf('%s J2 model in %.0f s\n', r1.label, toc);

        R(:,:,k) = gcrfRotation(cfg);
        x0(:,k) = 1e3 * [R(:,:,k)*out.r(:,1); R(:,:,k)*out.v(:,1)];
        altJ2(:,k) = orbitMean(vecnorm(out.r) - cfg.Re, out.T/dt);
        lowJ2(k) = min(vecnorm(out.r) - cfg.Re);
    end

    save(refFile, 'ref', 'R', 'x0', 'altJ2', 'lowJ2');
end

rsun = sunVector(t, ref(1).cfg); % rsun: 3 x N Sun position history, common to every case [km]
tt = datetime(ref(1).cfg.epoch) + seconds(t); % tt: 1 x N UTC output times

task = struct('label', {'EGM2008, no drag', 'EGM2008 + drag, F10.7 150', 'EGM2008 + drag, F10.7 150', ...
                        'EGM2008 + drag, CelesTrak', 'EGM2008 + drag, CelesTrak'}, ...
              'drag', {false, true, true, true, true}, ...
              'sw', {'none', 'constant', 'constant', 'file', 'file'}, ...
              'B', {1, 30, 150, 30, 150});
% task: toolbox configurations with drag flag, space weather source and assumed ballistic coefficient m/(Cd A) [kg/m^2]

swMat = ''; % swMat: space weather MAT file converted from swFile
if isfile(swFile)
    swMat = aeroReadSpaceWeatherData(swFile, FolderName=ckDir);
else
    task = task(~strcmp({task.sw}, 'file'));
    fprintf('%s not found, CelesTrak runs skipped\n', swFile);
end

nTask = numel(task); % nTask: number of toolbox configurations
nSeg = ceil(nDay/30); % nSeg: propagation segments per configuration, each saved as a checkpoint
edge = round(linspace(1, N, nSeg + 1)); % edge: 1 x nSeg+1 sample indices bounding the segments, shared at each boundary
P = cell(1, nTask); % P: 3 x N x nCase GCRF position histories of each configuration [m]
tRun = zeros(1, nTask); % tRun: wall time of each configuration in this session [s]

parfor j = 1:nTask
    opts = Aero.spacecraft.NumericalPropagatorOptions; % opts: numerical propagator options of this configuration
    opts.GravitationalPotentialModel = 'spherical-harmonics';
    opts.SphericalHarmonicDegree = deg;
    opts.ODESet = odeset('RelTol', tol, 'AbsTol', tol);
    opts.IncludeAtmosDrag = task(j).drag;
    if strcmp(task(j).sw, 'file')
        opts.UseSpaceWeatherDataFile = true;
        opts.SpaceWeatherDataFile = swMat;
    end

    pp = Aero.spacecraft.PhysicalProperties(Mass=task(j).B*ones(1, nCase), DragCoefficient=ones(1, nCase), ...
                                            DragArea=ones(1, nCase)); % pp: spacecraft with Cd = 1 and A = 1 m^2, so m/(Cd A) = B

    key = sprintf('%s_B%d_deg%d_tol%g_%dd', task(j).sw, task(j).B, deg, tol, nDay); % key: checkpoint name stem of this configuration
    Pj = zeros(3, N, nCase); % Pj: 3 x N x nCase GCRF position history of this configuration [m]
    r = x0(1:3,:); % r: 3 x nCase GCRF position at the segment start [m]
    v = x0(4:6,:); % v: 3 x nCase GCRF velocity at the segment start [m/s]

    t1 = tic; % t1: timer start of this configuration
    for s = 1:nSeg
        idx = edge(s) : edge(s+1); % idx: samples of this segment
        f = fullfile(ckDir, sprintf('%s_seg%02d.mat', key, s)); % f: checkpoint file of this segment

        if isfile(f)
            c = load(f); % c: segment positions p [m] and final velocity v [m/s]
        else
            [p, vs] = propagateOrbit(tt(idx), r, v, PropModel='numerical', ...
                                     NumericalPropagatorOptions=opts, PhysicalProperties=pp);
            % p: 3 x M x nCase segment positions [m], vs: 3 x M x nCase segment velocities [m/s]
            c = struct('p', p, 'v', reshape(vs(:,end,:), 3, nCase));
            saveSegment(f, c);
        end

        Pj(:, idx, :) = c.p;
        r = reshape(c.p(:,end,:), 3, nCase);
        v = c.v;
    end

    P{j} = Pj;
    tRun(j) = toc(t1);
end

summ = cell(1, nTask); % summ: plotTradeLighting.m summary of each configuration
alt = cell(1, nTask); % alt: N x nCase orbit averaged altitude of each configuration [km]
low = zeros(nTask, nCase); % low: lowest altitude of each configuration and case over the run [km]

f0 = findobj('Type', 'figure'); % f0: figures open before this script draws any
summJ2 = plotTradeLighting(ref); % summJ2: plotTradeLighting.m summary of the J2 model
saveFigs(outDir, sprintf('lighting_J2model_%dd', nDay), f0);

for j = 1:nTask
    res = ref; % res: this configuration in the tradeLighting.m result format
    alt{j} = zeros(N, nCase);

    for k = 1:nCase
        o.t = t; % o: state history struct in the propagate.m format
        o.r = R(:,:,k).' * P{j}(:,:,k) / 1e3;
        o.rsun = rsun;

        res(k).ecl = eclipseTimes(o, ref(k).cfg);
        alt{j}(:,k) = orbitMean(vecnorm(o.r) - ref(k).cfg.Re, ref(k).T/dt);
        low(j,k) = min(vecnorm(o.r) - ref(k).cfg.Re);
    end

    summ{j} = plotTradeLighting(res);
    saveFigs(outDir, sprintf('lighting_task%d_%dd', j, nDay), f0);
end

save(fullfile(outDir, sprintf('numericalCrossCheck_%dd.mat', nDay)), 'task', 'tRun', 'summ', 'summJ2', 'alt', 'altJ2', ...
     'low', 'lowJ2', 'deg', 'tol', 't');

fprintf('\ntoolbox wall time per configuration, %d cases together\n', nCase);
for j = 1:nTask
    fprintf('  %-28s B %3d   %6.0f s\n', task(j).label, task(j).B, tRun(j));
end

w = @(T) round(T/dt); % w: one orbital period in samples
fprintf('\ntoolbox minus J2 model over %d days, EGM2008 %d x %d, RelTol %.0e, B = m/(Cd A) assumed [kg/m^2]\n', nDay, deg, deg, tol);
fprintf('case  model full Sun [%%]  model longest [min]  model alt [km]   configuration               B    d full Sun [pp]  d longest [min]  d alt end [km]  decay [km]  lowest [km]\n');
for k = 1:nCase
    fprintf('%-4s  %18.2f  %19.2f  %14.2f   %-26s  %3s  %15s  %15s  %14s  %10s  %11.2f\n', ref(k).name, summJ2(k).fullSun, ...
            summJ2(k).longest, altJ2(w(ref(k).T), k), 'J2 model', '', '', '', '', '', lowJ2(k));
    for j = 1:nTask
        a1 = alt{j}(end - w(ref(k).T), k); % a1: orbit averaged toolbox altitude one orbit before the end [km]
        fprintf('%-4s  %18s  %19s  %14s   %-26s  %3d  %15.3f  %15.3f  %14.2f  %10.2f  %11.2f\n', ref(k).name, '', '', '', ...
                task(j).label, task(j).B*task(j).drag, summ{j}(k).fullSun - summJ2(k).fullSun, ...
                summ{j}(k).longest - summJ2(k).longest, a1 - altJ2(end - w(ref(k).T), k), a1 - alt{j}(w(ref(k).T), k), low(j,k));
    end
end

figure('Color', 'w', 'Position', [80 60 1100 750]);
for k = 1:nCase
    subplot(3, 2, k);
    plot(t/86400, altJ2(:,k) - altJ2(w(ref(k).T), k), 'k', 'LineWidth', 1.2); grid on; hold on;
    for j = 1:nTask
        plot(t/86400, alt{j}(:,k) - alt{j}(w(ref(k).T), k), 'LineWidth', 1.0);
    end
    xlim([0 nDay]);
    ylabel('\Delta orbit mean altitude [km]'); xlabel('days from epoch');
    title(ref(k).label, 'FontWeight', 'normal', 'FontSize', 9);
end
lg = [{'J2 model'}, arrayfun(@(s) sprintf('%s, B %d', s.label, s.B*s.drag), task, 'UniformOutput', false)]; % lg: legend entries
legend(lg, 'Location', 'southwest', 'FontSize', 7);
exportgraphics(gcf, fullfile(outDir, sprintf('numericalCrossCheck_altitude_%dd.png', nDay)), 'Resolution', 150);

diary off;


function y = orbitMean(x, n)
% orbitMean: running mean over one orbital period.
% x: 1 x N samples
% n: orbital period [samples]
% y: N x 1 running mean

y = movmean(x(:), round(n));

end


function saveFigs(outDir, stem, f0)
% saveFigs: saves the figures opened since f0 as PNG files and closes them.
% outDir: folder for figures
% stem: file name stem
% f0: figure handles to leave untouched

f = findobj('Type', 'figure'); % f: open figure handles
f = f(arrayfun(@(h) ~any(h == f0), f));
for j = 1:numel(f)
    exportgraphics(f(j), fullfile(outDir, sprintf('%s_%d.png', stem, f(j).Number)), 'Resolution', 120);
    close(f(j));
end

end


function saveSegment(f, c)
% saveSegment: writes one segment checkpoint under a temporary name, then renames it into place.
% f: checkpoint file
% c: struct of segment positions p [m] and final velocity v [m/s]

fTmp = strrep(f, '.mat', '_part.mat'); % fTmp: temporary file name used while writing
save(fTmp, '-struct', 'c');
movefile(fTmp, f, 'f');

end
