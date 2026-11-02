function fr = slewFrequency(out, cfg, D, tilt)
% slewFrequency: observation windows with the panels tilted away from the Sun, each started once the battery is back at the top of charge in Nominal, as many as the run allows, the panel temperature carried step by step.
% out: state history struct from propagate.m, cfg.use.eclipse on
% cfg: parameter struct from config.m
% D: 1 x L window lengths [min]
% tilt: 1 x M tilts of the panel normal from the Sun in a window, cfg.slew.tilt when omitted [deg]
% fr: M x L struct array, one per tilt and window length, of window starts and stops [s], recovery times [s], windows cut short by a drop out of Nominal, complete windows per day, lowest state of charge, time in each mode [s], and the soc, mode, tilt, array power and panel temperature histories

if nargin < 4
    tilt = cfg.slew.tilt;
end

N = numel(out.t); % N: number of output samples
days = (out.t(end) - out.t(1))/86400; % days: run length [days]

r = arrayPower(struct('t', out.t(1:2), 'illum', out.illum(1:2), 'rsun', out.rsun(:, 1:2)), cfg, cfg.power.panelArea, ...
               struct('incidenceAngles_deg', 0, 'incidenceSummaryAngles_deg', 0, 'warnTemperatureExtrapolation', false, 'warnLumpedModel', false));
% r: arrayPower over the first output step, run for its parameter struct
p = r.params; % p: arrayPower parameters, end-of-life factor from the sheet included

S = solarFlux(out, cfg); % S: 1 x N sunlight on a surface held normal to the Sun [W/m^2]
c0 = max(cosd(cfg.slew.tiltSun), 0); % c0: share of the sunlight on the front face outside the windows

for i = 1:numel(tilt)
    c1 = max(cosd(tilt(i)), 0); % c1: share of the sunlight on the front face in a window

    for j = 1:numel(D)
        soc = cfg.batt.socTop*ones(1, N); % soc: 1 x N state of charge, share of rated capacity
        mode = ones(1, N); % mode: 1 x N operating mode, 1 Nominal, 2 Raven, 3 Safe
        th = cfg.slew.tiltSun*ones(1, N); % th: 1 x N panel tilt from the Sun over the step from each sample [deg]
        pa = zeros(1, N); % pa: 1 x N array power at the bus, at the tilt of the step from each sample [W]
        T = (p.initialTemperature_C + 273.15)*ones(1, N); % T: 1 x N panel temperature, arrayPower's starting value first [K]
        t0 = zeros(1, 0); % t0: window starts [s]
        t1 = zeros(1, 0); % t1: window stops [s]
        cut = false(1, 0); % cut: windows stopped early by a drop out of Nominal
        on = false; % on: true while a window runs

        for k = 1:N-1
            if ~on && mode(k) == 1 && soc(k) >= cfg.batt.socTop
                on = true;
                t0(end+1) = out.t(k);
            end
            if on && (out.t(k) >= t0(end) + 60*D(j) || mode(k) ~= 1)
                on = false;
                t1(end+1) = out.t(k);
                cut(end+1) = mode(k) ~= 1;
            end

            cs = c0; % cs: share of the sunlight on the front face over the step
            if on
                th(k) = tilt(i);
                cs = c1;
            end

            [T(k+1), pS, pE] = panelStep(T(k), S(k:k+1)*cs, out.t(k+1) - out.t(k), p, cfg.power.panelArea);
            % pS: array power at the bus at the start of the step [W]
            % pE: array power at the bus at the end of the step [W]
            pa(k) = pS;
            pa(k+1) = pE;

            % Array power at its mean over the step:
            [soc(k+1), mode(k+1)] = modeStep(soc(k), mode(k), (pS + pE)/2, out.t(k+1) - out.t(k), cfg);
        end
        if on
            t1(end+1) = out.t(N);
            cut(end+1) = false;
        end

        rec = NaN(1, numel(t1)); % rec: time from each window stop back to the top of charge, NaN if not within the run [s]
        for q = 1:numel(t1)
            kq = find(out.t >= t1(q) & soc >= cfg.batt.socTop, 1); % kq: first sample at the top of charge after the stop
            if ~isempty(kq)
                rec(q) = out.t(kq) - t1(q);
            end
        end

        done = ~cut & (t1 - t0 >= 60*D(j)); % done: windows run for their full length

        fr(i, j).tilt = tilt(i); % fr.tilt: panel tilt from the Sun in the windows [deg]
        fr(i, j).D = D(j); % fr.D: window length [min]
        fr(i, j).t0 = t0; % fr.t0: window starts [s]
        fr(i, j).t1 = t1; % fr.t1: window stops [s]
        fr(i, j).cut = cut; % fr.cut: windows stopped early by a drop out of Nominal
        fr(i, j).rec = rec; % fr.rec: time from each window stop back to the top of charge [s]
        fr(i, j).perDay = nnz(done)/days; % fr.perDay: windows run for their full length per day
        fr(i, j).socMin = min(soc); % fr.socMin: lowest state of charge, share of rated capacity
        fr(i, j).tMode = accumarray(mode(1:N-1).', diff(out.t).', [3 1]).'; % fr.tMode: 1 x 3 time in Nominal, Raven and Safe [s]
        fr(i, j).soc = soc; % fr.soc: 1 x N state of charge
        fr(i, j).mode = mode; % fr.mode: 1 x N operating mode
        fr(i, j).th = th; % fr.th: 1 x N panel tilt from the Sun over the step from each sample [deg]
        fr(i, j).pArr = pa; % fr.pArr: 1 x N array power at the bus [W]
        fr(i, j).T = T; % fr.T: 1 x N panel temperature [K]
    end
end

end
