% numericalSettings: one day runtime and accuracy of toolbox numerical propagator settings for the six trade cases.

here = fileparts(mfilename('fullpath')); % here: folder holding this script
addpath(fullfile(here, '..'));
odinPath;

outDir = fullfile(here, 'results'); % outDir: folder for logs and figures
if ~isfolder(outDir)
    mkdir(outDir);
end
logFile = fullfile(outDir, 'numericalSettings.txt'); % logFile: text log of this run
if isfile(logFile)
    delete(logFile);
end
diary(logFile);

nDay = 1; % nDay: run length [days]
nCase = 6; % nCase: number of trade cases

x0 = zeros(6, nCase); % x0: initial GCRF state of each case [m; m/s]
rJ2 = []; % rJ2: 3 x N x nCase GCRF position history of the J2 model [m]
lbl = cell(nCase, 1); % lbl: case labels

for k = 1:nCase
    res = tradeLighting(nDay, k); % res: case parameters from tradeLighting.m
    cfg = res.cfg; % cfg: parameter struct of this case
    lbl{k} = res.label;

    out = propagate(cfg, nDay*86400/res.T); % out: J2 model state history on the 60 s grid
    R = gcrfRotation(cfg); % R: mean equator and equinox of epoch to GCRF rotation

    x0(:,k) = 1e3 * [R*out.r(:,1); R*out.v(:,1)];
    rJ2(:,:,k) = 1e3 * R*out.r;
end

tt = datetime(cfg.epoch) + seconds(out.t); % tt: 1 x N UTC output times

setting = {'J2 only, RelTol 1e-8',           'oblate-ellipsoid',    120, 1e-8,  false
           'EGM2008 120 x 120, RelTol 1e-10', 'spherical-harmonics', 120, 1e-10, false
           'EGM2008 120 x 120, RelTol 1e-8',  'spherical-harmonics', 120, 1e-8,  false
           'EGM2008 36 x 36, RelTol 1e-8',    'spherical-harmonics',  36, 1e-8,  false
           'EGM2008 20 x 20, RelTol 1e-8',    'spherical-harmonics',  20, 1e-8,  false
           'EGM2008 120 x 120 + drag, 1e-8',  'spherical-harmonics', 120, 1e-8,  true};
% setting: label, gravity model, degree and order, relative tolerance and drag flag of each setting

nRun = size(setting, 1); % nRun: number of settings
P = zeros([size(rJ2), nRun]); % P: 3 x N x nCase x nRun GCRF position histories [m]
tRun = zeros(nRun, 1); % tRun: wall time of each setting [s]

pp = Aero.spacecraft.PhysicalProperties(Mass=30*ones(1, nCase), DragCoefficient=ones(1, nCase), DragArea=ones(1, nCase));
% pp: spacecraft with ballistic coefficient m/(Cd A) = 30 kg/m^2

for j = 1:nRun
    opts = Aero.spacecraft.NumericalPropagatorOptions; % opts: numerical propagator options of this setting
    opts.GravitationalPotentialModel = setting{j,2};
    opts.ODESet = odeset('RelTol', setting{j,4}, 'AbsTol', setting{j,4});
    opts.IncludeAtmosDrag = setting{j,5};
    if strcmp(setting{j,2}, 'spherical-harmonics')
        opts.SphericalHarmonicDegree = setting{j,3};
    end

    tic;
    P(:,:,:,j) = propagateOrbit(tt, x0(1:3,:), x0(4:6,:), PropModel='numerical', ...
                                NumericalPropagatorOptions=opts, PhysicalProperties=pp);
    tRun(j) = toc;
end

fprintf('\nwall time for %d day, %d cases propagated together\n', nDay, nCase);
for j = 1:nRun
    fprintf('  %-34s %7.1f s, %6.1f h for a year\n', setting{j,1}, tRun(j), tRun(j)*365/nDay/3600);
end

fprintf('\nlargest position difference after %d day [m]\n', nDay);
fprintf('case   toolbox J2 - model J2   120 at 1e-8 - 1e-10   36 - 120 at 1e-8   20 - 120 at 1e-8   drag - no drag\n');
for k = 1:nCase
    fprintf('%-4s   %21.2f   %19.2f   %16.2f   %16.2f   %14.2f\n', lbl{k}(1), ...
            max(vecnorm(P(:,:,k,1) - rJ2(:,:,k))), max(vecnorm(P(:,:,k,3) - P(:,:,k,2))), ...
            max(vecnorm(P(:,:,k,4) - P(:,:,k,3))), max(vecnorm(P(:,:,k,5) - P(:,:,k,3))), ...
            max(vecnorm(P(:,:,k,6) - P(:,:,k,3))));
end

diary off;
save(fullfile(outDir, 'numericalSettings.mat'), 'setting', 'tRun', 'P', 'rJ2', 'x0', 'lbl');
