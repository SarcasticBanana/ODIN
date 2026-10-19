% groundContacts: ground station contacts, a first UHF link budget and contact time inside shadow for the six trade cases.

here = fileparts(which('groundContacts')); % here: folder holding this script
addpath(fullfile(here, '..'));
odinPath;

outDir = fullfile(here, 'results'); % outDir: folder for logs, tables and figures
if ~isfolder(outDir)
    mkdir(outDir);
end
logFile = fullfile(outDir, 'groundContacts.txt'); % logFile: text log of this run
if isfile(logFile)
    delete(logFile);
end
diary(logFile);

gsName = 'Boulder campus placeholder'; % gsName: station label, a stand-in until Comms picks sites
gsLat = 40.007; % gsLat: station geodetic latitude [deg]
gsLon = -105.266; % gsLon: station longitude, east positive [deg]
gsAlt = 0; % gsAlt: station height above the ellipsoid [m]
elMin = 10; % elMin: elevation mask [deg]
f = 437e6; % f: downlink frequency [Hz]
nDay = 365; % nDay: mission duration [days]
dt = 60; % dt: output and scenario sample spacing [s]

Ptx = 0; % Ptx: placeholder transmit power, 1 W [dBW]
Gtx = 0; % Gtx: placeholder spacecraft antenna gain toward the station [dBi]
Ltx = 1; % Ltx: placeholder spacecraft line loss [dB]
GT = -13; % GT: placeholder station G/T, 14 dBi Yagi over 500 K [dB/K]
Lpol = 3; % Lpol: placeholder polarisation mismatch loss [dB]
Latm = 1; % Latm: placeholder atmospheric and ionospheric allowance [dB]
Rb = 9600; % Rb: placeholder data rate [bit/s]
EbN0req = 10; % EbN0req: placeholder required Eb/N0 [dB]
c0 = 299792458; % c0: speed of light [m/s]
kB = -228.6; % kB: Boltzmann constant [dBW/(K Hz)]

% Link margin with free-space loss at slant range r, Eb/N0 = EIRP - Lfs - L + G/T - k - 10 log10(Rb):
margin = @(r) Ptx + Gtx - Ltx - 20*log10(4*pi*r*f/c0) - Lpol - Latm + GT - kB - 10*log10(Rb) - EbN0req;
% margin: link margin at a slant range given in metres [dB]

nCase = 6; % nCase: number of trade cases
S = struct('name', {}, 'nPerDay', {}, 'durMean', {}, 'durMax', {}, 'minPerDay', {}, 'mMask', {}, 'mBest', {}, ...
           'volPerDay', {}, 'nShadow', {}, 'shadowMinPerDay', {}, 'shadowMax', {}, 'umbraMax', {}); % S: per case summary

fg = figure('Color', 'w', 'Position', [80 60 1100 750], 'Visible', 'off'); % fg: contact summary figure, drawn off screen

ckDir = fullfile(tempdir, 'odin_groundContacts'); % ckDir: checkpoint folder, so an interrupted run resumes where it stopped
if ~isfolder(ckDir)
    mkdir(ckDir);
end

for k = 1:nCase
    ckFile = fullfile(ckDir, sprintf('case%d_%dd.mat', k, nDay)); % ckFile: checkpoint of this case
    if isfile(ckFile)
        load(ckFile, 'Sk', 'aos', 'dur', 'inShadow', 'lbl');
        fprintf('%s loaded from %s\n', lbl, ckFile);
    else
        r1 = tradeLighting(1/24, k); % r1: case parameters from tradeLighting.m, one hour run
        cfg = r1.cfg; % cfg: parameter struct of this case

        out = propagate(cfg, (nDay*86400 + 1)/r1.T); % out: model state history on the 60 s grid, one second past the last output time
        ecl = eclipseTimes(out, cfg); % ecl: exact shadow entry and exit times [s]
        [p0, p1] = closeIntervals(ecl.penumbraIn, ecl.penumbraOut, out.t(end)); % p0, p1: shadow interval starts and ends [s]
        [u0, u1] = closeIntervals(ecl.umbraIn, ecl.umbraOut, out.t(end)); % u0, u1: umbra interval starts and ends [s]

        t0 = datetime(cfg.epoch); % t0: UTC epoch
        sc = satelliteScenario(t0, t0 + seconds(nDay*86400), dt); % sc: scenario sampled on the output grid
        sat = odinToScenario(sc, out, r1.name); % sat: ODIN ephemeris satellite of this case
        gs = groundStation(sc, gsLat, gsLon, 'Altitude', gsAlt, 'MinElevationAngle', elMin, 'Name', gsName); % gs: ground station

        iv = accessIntervals(access(sat, gs)); % iv: toolbox contact intervals, edges on the sample grid
        [~, elH, rH] = aer(gs, sat); % elH: 1 x N elevation [deg], rH: 1 x N slant range [m], both at the sample times
        tH = (0:numel(elH)-1) * dt; % tH: 1 x N sample times since epoch [s]

        i0 = round(seconds(iv.StartTime - sc.StartTime).'/dt) + 1; % i0: 1 x M first sample of each contact
        i1 = round(seconds(iv.EndTime - sc.StartTime).'/dt) + 1; % i1: 1 x M last sample of each contact
        M = numel(i0); % M: number of contacts

        aos = tH(i0); % aos: 1 x M refined acquisition times [s]
        los = tH(i1); % los: 1 x M refined loss of signal times [s]
        elMax = zeros(1, M); % elMax: 1 x M highest elevation of each contact [deg]
        mMask = zeros(1, M); % mMask: 1 x M link margin at the elevation mask [dB]
        mBest = zeros(1, M); % mBest: 1 x M link margin at the highest elevation [dB]
        vol = zeros(1, M); % vol: 1 x M data volume with positive margin [bit]
        inShadow = zeros(1, M); % inShadow: 1 x M contact time inside any shadow [s]
        inUmbra = zeros(1, M); % inUmbra: 1 x M contact time inside the umbra [s]

        for m = 1:M
            w = max(1, i0(m) - 3) : min(numel(tH), i1(m) + 3); % w: samples spanning the contact and three either side
            ppEl = spline(tH(w), elH(w)); % ppEl: elevation spline across the contact [deg]
            ppR = spline(tH(w), rH(w)); % ppR: slant range spline across the contact [m]

            if i0(m) > 1
                aos(m) = fzero(@(x) ppval(ppEl, x) - elMin, [tH(i0(m) - 1), tH(i0(m))]);
            end
            if i1(m) < numel(tH)
                los(m) = fzero(@(x) ppval(ppEl, x) - elMin, [tH(i1(m)), tH(i1(m) + 1)]);
            end

            [xMax, negEl] = fminbnd(@(x) -ppval(ppEl, x), aos(m), los(m)); % xMax: time of the highest elevation [s], negEl: minus that elevation [deg]
            elMax(m) = -negEl;
            mMask(m) = margin(ppval(ppR, aos(m)));
            mBest(m) = margin(ppval(ppR, xMax));

            tq = linspace(aos(m), los(m), max(2, ceil((los(m) - aos(m))/10))); % tq: sample times across the contact [s]
            vol(m) = Rb * trapz(tq, double(margin(ppval(ppR, tq)) >= 0));

            % Interval overlap, sum of max(0, min(ends) - max(starts)):
            inShadow(m) = sum(max(0, min(p1, los(m)) - max(p0, aos(m))));
            inUmbra(m) = sum(max(0, min(u1, los(m)) - max(u0, aos(m))));
        end

        dur = los - aos; % dur: 1 x M contact durations [s]
        Sk = struct('name', r1.name, 'nPerDay', M/nDay, 'durMean', mean(dur)/60, 'durMax', max(dur)/60, ...
                    'minPerDay', sum(dur)/60/nDay, 'mMask', min(mMask), 'mBest', max(mBest), 'volPerDay', sum(vol)/8/1e6/nDay, ...
                    'nShadow', nnz(inShadow > 0), 'shadowMinPerDay', sum(inShadow)/60/nDay, 'shadowMax', max([0 inShadow])/60, ...
                    'umbraMax', max([0 inUmbra])/60); % Sk: summary of this case

        T = table(string(sc.StartTime + seconds(aos.'), 'yyyy-MM-dd HH:mm:ss.SSS'), string(sc.StartTime + seconds(los.'), 'yyyy-MM-dd HH:mm:ss.SSS'), ...
                  dur.'/60, elMax.', mMask.', mBest.', vol.'/8/1e3, inShadow.'/60, inUmbra.'/60, ...
                  'VariableNames', {'AOS_UTC', 'LOS_UTC', 'duration_min', 'maxElevation_deg', 'marginAtMask_dB', 'marginAtMaxEl_dB', ...
                                    'volume_kB', 'inShadow_min', 'inUmbra_min'}); % T: contact table of this case
        writetable(T, fullfile(outDir, sprintf('groundContacts_case%s.csv', r1.name)));

        lbl = r1.label; % lbl: case label
        save(ckFile, 'Sk', 'aos', 'dur', 'inShadow', 'lbl');
        fprintf('%s: %d contacts\n', lbl, M);
    end

    S(k) = Sk;

    ax = subplot(3, 2, k, 'Parent', fg); % ax: axes of this case
    day = floor(aos/86400); % day: 1 x M day of each contact since epoch
    plot(ax, 0:nDay-1, accumarray(day(:)+1, dur(:)/60, [nDay 1]), 'LineWidth', 1.0); grid(ax, 'on'); hold(ax, 'on');
    plot(ax, 0:nDay-1, accumarray(day(:)+1, inShadow(:)/60, [nDay 1]), 'LineWidth', 1.0);
    xlim(ax, [0 nDay]);
    ylabel(ax, 'contact per day [min]'); xlabel(ax, 'days from epoch');
    title(ax, lbl, 'FontWeight', 'normal', 'FontSize', 9);
    if k == 1
        legend(ax, 'all contacts', 'inside shadow', 'Location', 'best');
    end
end

exportgraphics(fg, fullfile(outDir, 'groundContacts.png'), 'Resolution', 150);
close(fg);

fprintf('\n%s, %.3f N %.3f E, %g deg mask, %g MHz, %d days\n', gsName, gsLat, gsLon, elMin, f/1e6, nDay);
fprintf('placeholder radio: %g dBW, %g dBi, %g dB line, G/T %g dB/K, %g dB polarisation, %g dB atmosphere, %g bit/s, Eb/N0 %g dB\n', ...
        Ptx, Gtx, Ltx, GT, Lpol, Latm, Rb, EbN0req);
fprintf('case  contacts/day  mean [min]  longest [min]  min/day  margin at mask [dB]  best margin [dB]  MB/day   in shadow: contacts  min/day  longest [min]  longest in umbra [min]\n');
for k = 1:nCase
    fprintf('%-4s  %12.2f  %10.2f  %13.2f  %7.1f  %19.1f  %16.1f  %6.2f   %18d  %7.2f  %13.2f  %22.2f\n', S(k).name, S(k).nPerDay, ...
            S(k).durMean, S(k).durMax, S(k).minPerDay, S(k).mMask, S(k).mBest, S(k).volPerDay, S(k).nShadow, ...
            S(k).shadowMinPerDay, S(k).shadowMax, S(k).umbraMax);
end

diary off;
save(fullfile(outDir, 'groundContacts.mat'), 'S', 'gsName', 'gsLat', 'gsLon', 'gsAlt', 'elMin', 'f', 'nDay');


function [t0, t1] = closeIntervals(tIn, tOut, tEnd)
% closeIntervals: entry and exit times paired into intervals, closed at the run start and end.
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
