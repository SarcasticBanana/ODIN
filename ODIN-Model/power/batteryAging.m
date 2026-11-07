function b = batteryAging(cfg)
% batteryAging: battery capacity fade and resistance growth over the mission, one eclipse every orbit, NMC and LFP cell models of NREL BLAST-Lite.
% cfg: parameter struct from config.m, battery inputs in cfg.batt and the mission length in cfg.power.years
% b: struct of daily histories t [day], efc, q, r, dod and socMin, and the values eDraw [Wh], dodBol, kR, qEOL, rEOL, socMinEOL, etaDisEOL, yearsFloor [yr] and tFail [day]

c = cfg.batt; % c: battery inputs, Battery section of the Inputs sheet

% Kepler's third law:
tOrb = 2*pi*sqrt((cfg.Re + c.h)^3/cfg.mu)/60; % tOrb: period of the battery sizing orbit [min]

% Battery energy per eclipse, bus load through the longest shadow over the discharge efficiency:
b.eDraw = c.pEcl*c.tEcl/60/c.etaDis; % b.eDraw: battery energy per eclipse at the start of life [Wh]
b.dodBol = b.eDraw/c.energy; % b.dodBol: depth of discharge per eclipse at the start of life, share of rated capacity

% Battery current, I = E/(t V):
I = b.eDraw/(c.tEcl/60)/c.vNom; % I: battery current through the shadow [A]

% Pack resistance, R = ns rho/(E/V), with rho the cell resistance times the cell capacity:
R0 = c.ns*c.rho/(c.energy/c.vNom); % R0: pack resistance at the start of life, cells only [Ohm]

% I^2 R loss over the shadow per unit of resistance growth, as a share of the battery energy per eclipse:
b.kR = I^2*R0*(c.tEcl/60)/b.eDraw; % b.kR: deepening of the eclipse per unit of resistance growth

s = agingRun(c.model, b.dodBol, c.socTop, c.temp, cfg.power.years, tOrb, c.tEcl, c.fadeScale, c.dose, b.kR);
% s: daily histories from agingRun

b.t = s.t; % b.t: 1 x N days from the start of the mission [day]
b.efc = s.efc; % b.efc: 1 x N equivalent full cycles so far
b.q = s.q; % b.q: 1 x N capacity over its start-of-life value
b.r = s.r; % b.r: 1 x N resistance over its start-of-life value
b.dod = s.dod; % b.dod: 1 x N depth of discharge per eclipse, share of the remaining capacity

% Lowest state of charge, share of full charge (the energy held at the top of charge), 1 - E_draw (1 + kR (r - 1))/(q E socTop):
b.socMin = 1 - b.eDraw*(1 + b.kR*(b.r - 1))./(b.q*c.energy*c.socTop); % b.socMin: 1 x N lowest state of charge in the eclipse, share of full charge

b.qEOL = b.q(end); % b.qEOL: capacity at the end of the mission over its start-of-life value
b.rEOL = b.r(end); % b.rEOL: resistance at the end of the mission over its start-of-life value
b.tFail = s.tFail; % b.tFail: first day the eclipse needs more charge than the battery holds, NaN if never [day]
b.socMinEOL = b.socMin(end); % b.socMinEOL: lowest state of charge at the end of the mission, NaN if the battery runs empty
if ~isnan(b.tFail)
    b.socMinEOL = NaN;
end

% Discharge efficiency with the grown resistance, eta/(1 + kR (r - 1)):
b.etaDisEOL = c.etaDis/(1 + b.kR*(b.rEOL - 1)); % b.etaDisEOL: discharge efficiency at the end of the mission

k = find(b.socMin < c.floor - 1e-12, 1); % k: first day below the lowest allowed state of charge
b.yearsFloor = NaN; % b.yearsFloor: years until the lowest state of charge drops below cfg.batt.floor, NaN if never [yr]
if ~isempty(k)
    b.yearsFloor = b.t(k)/365.25;
end

end


function s = agingRun(model, dod0, socTop, TC, years, tOrb, tEcl, scale, dose, kR)
% agingRun: daily capacity and resistance with the same battery energy every eclipse, rates re-read each day at the current depth of discharge.
% model: cell aging model, 'NMC' or 'LFP'
% dod0: depth of discharge per eclipse at the start of life, share of rated capacity
% socTop: top of charge, share of rated capacity
% TC: cell temperature [deg C]
% years: mission length [yr]
% tOrb: orbital period [min]
% tEcl: eclipse length [min]
% scale: factor on every capacity-loss and resistance-growth term, 1 for the published model
% dose: total ionising dose at the cells over the mission [krad]
% kR: deepening of the eclipse per unit of resistance growth
% s: struct of 1 x N daily histories t [day], efc, q, r and dod, and tFail [day], NaN if the battery never runs empty

if dod0 >= socTop
    s = struct('t', 0, 'efc', 0, 'q', 1, 'r', 1, 'dod', dod0, 'tFail', 0);
    return
end

nDay = round(years*365.25); % nDay: days in the mission
nOrb = 1440/tOrb; % nOrb: eclipses per day

% Rate table over depth of discharge:
g = linspace(max(0.01, 0.5*dod0), min(socTop, 1), 120)'; % g: 120 x 1 depths of discharge of the table, share of rated capacity
cr = 2*dod0/(tOrb/60); % cr: orbit-average C-rate on rated capacity, fixed by the fixed charge per eclipse [1/h]
R = orbitRates(model, g, socTop, TC, tOrb, tEcl, cr); % R: struct of 120 x 1 rates from orbitRates

s.t = 0:nDay; % s.t: 1 x N days [day]
s.efc = zeros(1, nDay + 1); % s.efc: equivalent full cycles so far
s.q = ones(1, nDay + 1); % s.q: capacity over its start-of-life value
s.r = ones(1, nDay + 1); % s.r: resistance over its start-of-life value
s.dod = dod0*ones(1, nDay + 1); % s.dod: depth of discharge per eclipse, share of the remaining capacity
s.tFail = NaN; % s.tFail: first day the eclipse needs more charge than the battery holds [day]

qc = 0; % qc: calendar capacity loss before scaling
qy = 0; % qy: cycling capacity loss before scaling
qb = 0; % qb: LFP break-in capacity loss before scaling
rc = 0; % rc: calendar resistance growth before scaling
ry = 0; % ry: cycling resistance growth before scaling
n = nDay + 1; % n: length of the histories

for day = 1:nDay
    d = min(max(s.dod(day), g(1)), g(end)); % d: depth of discharge at the start of the day, held inside the table
    at = @(x) interp1(g, x, d); % at: rate at d from its column of the table

    % Equivalent full cycles in a day, each eclipse deepened by the I^2 R loss:
    dEfc = nOrb*dod0*(1 + kR*(s.r(day) - 1)); % dEfc: equivalent full cycles in the day

    switch model
        case 'NMC'
            % Charge throughput, two half cycles of the 2.15 Ah Sanyo cell per equivalent full cycle:
            ah = 2*2.15*dEfc; % ah: charge throughput in the day [Ah]

            % Schmalstieg: calendar alpha t^0.75 with t in days, cycling beta Ah^0.5, resistance alpha t^0.75 + beta Ah:
            qc = qc + dPower(qc, 1, at(R.aCap), 0.75);
            qy = qy + dPower(qy, ah, at(R.bCap), 0.5);
            rc = rc + dPower(rc, 1, at(R.aRes), 0.75);
            ry = ry + at(R.bRes)*ah;
        case 'LFP'
            % Naumann: calendar sigmoid, cycling (k EFC)^p, break-in sigmoid above 3 cycles a day, linear resistance:
            qc = qc + dSigmoid(qc, 1, at(R.q1), 1.30510034211874e-4, at(R.q3));
            qy = qy + dPowerB(qy, dEfc, at(R.q5), 1.12847759334355);
            if dEfc > 3
                qb = qb + dSigmoid(qb, dEfc, at(R.q7), 3.03553871631028e-3, 1.43752162947637);
            end
            rc = rc + at(R.rCal)*86400;
            ry = ry + at(R.rCyc)*dEfc/100;
    end

    % Radiation, 8.4% capacity per 2744 krad (Qiu et al., cited by Ali et al. 2023), linear over the mission:
    qRad = 0.084*dose/2744*(day/nDay); % qRad: capacity lost to radiation so far

    s.efc(day + 1) = s.efc(day) + dEfc;
    s.q(day + 1) = 1 - scale*(qc + qy + qb) - qRad;
    s.r(day + 1) = 1 + scale*(rc + ry);

    % Same battery energy every eclipse, dod = dod0 (1 + kR (r - 1))/q:
    s.dod(day + 1) = dod0*(1 + kR*(s.r(day + 1) - 1))/max(s.q(day + 1), 1e-6);

    if s.dod(day + 1) >= socTop
        s.tFail = day;
        n = day + 1;
        break
    end
end

s.t = s.t(1:n);
s.efc = s.efc(1:n);
s.q = s.q(1:n);
s.r = s.r(1:n);
s.dod = s.dod(1:n);

end


function R = orbitRates(model, dod, socTop, TC, tOrb, tEcl, cr)
% orbitRates: aging rates averaged over one orbit's state of charge profile, as BLAST-Lite does, at each depth of discharge.
% model: cell aging model, 'NMC' (Schmalstieg et al. 2014) or 'LFP' (Naumann et al. 2018 and 2020, NREL identification)
% dod: n x 1 depths of discharge per eclipse, share of rated capacity
% socTop: top of charge, share of rated capacity
% TC: cell temperature [deg C]
% tOrb: orbital period [min]
% tEcl: eclipse length [min]
% cr: orbit-average C-rate on rated capacity [1/h]
% R: struct of n x 1 rates, aCap, bCap, aRes and bRes for NMC, q1, q3, q5, q7, rCal and rCyc for LFP

T = tOrb*60; % T: orbital period [s]
te = tEcl*60; % te: eclipse length [s]
tc = 0.6*(T - te); % tc: recharge time, 60% of the sunlit time [s]
t = 0:10:(T + 10); % t: 1 x M sample times 10 s apart, the last one the first at or after T [s]
t = t(t < T + 10);
TK = TC + 273.15; % TK: cell temperature [K]

kE = t <= te; % kE: samples in the shadow
kC = t > te & t <= te + tc; % kC: samples recharging

% Linear discharge through the shadow, linear recharge over tc, then held at the top of charge:
soc = socTop*ones(numel(dod), numel(t)); % soc: n x M state of charge over one orbit, share of rated capacity
soc(:, kE) = socTop - dod*t(kE)/te;
soc(:, kC) = socTop - dod + dod*(t(kC) - te)/tc;

d = max(soc, [], 2) - min(soc, [], 2); % d: n x 1 depth of discharge of the sampled profile

% Time average over the orbit by the trapezoidal rule:
avg = @(y) trapz(t, y, 2)/(t(end) - t(1)); % avg: n x 1 orbit averages of the rows of y

switch model
    case 'NMC'
        V = nmcOcv(soc); % V: n x M open-circuit voltage of the Sanyo UR18650E cell [V]
        Vrms = sqrt(mean(V.^2, 2)); % Vrms: n x 1 root mean square of V over the samples [V]

        % Schmalstieg calendar capacity rate, (7.543 V - 23.75)e6 exp(-6976/T):
        R.aCap = avg((7.543*V - 23.75)*1e6*exp(-6976/TK));

        % Schmalstieg cycling capacity rate, 7.348e-3 (Vrms - 3.667)^2 + 7.6e-4 + 4.081e-3 DoD:
        R.bCap = 7.348e-3*(Vrms - 3.667).^2 + 7.6e-4 + 4.081e-3*d;

        % Schmalstieg calendar resistance rate, (5.270 V - 16.32)e5 exp(-5986/T):
        R.aRes = avg((5.270*V - 16.32)*1e5*exp(-5986/TK));

        % Schmalstieg cycling resistance rate, 2.153e-4 (Vrms - 3.725)^2 - 1.521e-5 + 2.798e-4 DoD:
        R.bRes = 2.153e-4*(Vrms - 3.725).^2 - 1.521e-5 + 2.798e-4*d;
    case 'LFP'
        Ua = uaGraphite(soc); % Ua: n x M graphite anode potential against Li/Li+ [V]

        % Calendar fade limit and time exponent, NREL fit to Naumann et al. 2018:
        R.q1 = avg(abs(0.989687151293590*exp(-2881067.56019324/TK^2*Ua.^0.5).*exp(8742.06309157261/TK*Ua.^0.5)));
        R.q3 = avg(abs(3.32850281062177e-4*exp(734553185711.369/TK^4*Ua.^(1/3)).*exp(-2.82161575620780e-6*TK^3*Ua.^0.25) ...
            .*exp(-3284991315.45121/TK^3*Ua.^(1/3)).*exp(1.27227593657290e-3*TK^2*Ua.^0.25)));

        % Cycling fade rate, NREL fit to Naumann et al. 2020:
        R.q5 = abs(-6.81260579372875e-6 + 2.59615973160844e-5*d + 2.11559710307295e-6*exp(d.^2*cr^3));

        % Break-in fade limit, skew-normal in state of charge and depth of discharge times a sigmoid in depth of discharge:
        R.q7 = avg(abs(0.582258029148225*skewNormal(soc, 0.0583128906965484, 0.208738181522897) ...
            .*skewNormal(d, -3.80744333129564, 1.16126260428210).*(2*(0.5 - 1./(1 + exp(25.4130804598602*d))))));

        % Calendar resistance rate, Arrhenius in temperature, cubic in state of charge:
        R.rCal = avg(3.4194e-10*exp(-71827/8.3144*(1/TK - 1/298.15))*(-3.3903*(soc - 0.5).^3 + 1.5604));

        % Cycling resistance rate, linear in C-rate, cubic in depth of discharge:
        R.rCyc = (-0.002*cr + 0.0021)*(6.8477*(d - 0.5).^3 + 0.91882);
    otherwise
        error('batteryAging:model', 'cfg.batt.model must be ''NMC'' or ''LFP'', not ''%s''', model);
end

end


function V = nmcOcv(soc)
% nmcOcv: open-circuit voltage of the Sanyo UR18650E cell, table from BLAST-Lite (Ecker et al. 2014, Schmalstieg et al. 2014).
% soc: state of charge, share of rated capacity, any size
% V: open-circuit voltage, same size as soc [V]

x = [0 0.008637153 0.026779514 0.044921875 0.063064236 0.081206597 0.099348958 0.117491319 ...
     0.135633681 0.153776042 0.171918403 0.190060764 0.208203125 0.226345486 0.244487847 0.262630208 ...
     0.280772569 0.298914931 0.317057292 0.335199653 0.353342014 0.371484375 0.389626736 0.407769097 ...
     0.425911458 0.444053819 0.462196181 0.480338542 0.498480903 0.516623264 0.534765625 0.552907986 ...
     0.571050347 0.589192708 0.607335069 0.625477431 0.643619792 0.661762153 0.679904514 0.698046875 ...
     0.716189236 0.734331597 0.752473958 0.770616319 0.788758681 0.806901042 0.825043403 0.843185764 ...
     0.861328125 0.879470486 0.897612847 0.915755208 0.933897569 0.952039931 0.970182292 0.988324653 ...
     0.998220486 1]; % x: state of charge of the table
y = [3.331 3.345014187 3.37917149 3.411603677 3.440585632 3.466289865 3.490268982 3.511315401 ...
     3.529946658 3.547197821 3.561688798 3.574972194 3.586357962 3.597053683 3.605506753 3.613442288 ...
     3.620342753 3.626380661 3.632073544 3.63690387 3.642079219 3.646909545 3.652084894 3.657605266 ...
     3.663470662 3.670198615 3.677444104 3.686759732 3.696420384 3.708323686 3.720572012 3.734372943 ...
     3.749553967 3.765252525 3.781123596 3.797857224 3.814590852 3.832532062 3.849783226 3.866861877 ...
     3.884458064 3.900846669 3.917752809 3.934831461 3.953462717 3.971921462 3.991415276 4.011771649 ...
     4.03195551 4.05196686 4.070770628 4.087849279 4.104237884 4.120108955 4.135980025 4.152541142 ...
     4.160649189 4.162]; % y: open-circuit voltage of the table [V]

V = interp1(x, y, min(max(soc, 0), 1));

end


function Ua = uaGraphite(soc)
% uaGraphite: graphite anode potential against Li/Li+ from the state of charge (Safari and Delacourt 2011, as in BLAST-Lite).
% soc: state of charge, share of rated capacity, any size
% Ua: anode potential, same size as soc [V]

x = 8.5e-3 + soc*(0.78 - 8.5e-3); % x: lithium fraction of the graphite

Ua = 0.6379 + 0.5416*exp(-305.5309*x) + 0.044*tanh(-(x - 0.1958)/0.1088) - 0.1978*tanh((x - 1.0571)/0.0854) ...
    - 0.6875*tanh((x + 0.0117)/0.0529) - 0.0175*tanh((x - 0.5692)/0.0875);

end


function f = skewNormal(x, a, w)
% skewNormal: skew-normal density centred on 0.5, as in the BLAST-Lite LFP break-in term.
% x: argument, any size
% a: skewness
% w: width
% f: density, same size as x

z = (x - 0.5)/w; % z: standardised argument

% 2 phi(z) Phi(a z):
f = 2*(exp(-z.^2/2)/sqrt(2*pi))*0.5.*(1 + erf(a*z/sqrt(2)));

end


function dy = dPower(y, dx, k, p)
% dPower: change of y = k x^p over a step dx from the x that gives the current y, the BLAST-Lite state update.
% y: current value
% dx: step in x
% k: rate
% p: exponent
% dy: change in y

if dx == 0
    dy = 0;
elseif y <= 0
    dy = k*dx^p;
else
    % Slope of k x^p at x = (y/k)^(1/p):
    dy = k*p*(y/k)^((p - 1)/p)*dx;
end

end


function dy = dPowerB(y, dx, k, p)
% dPowerB: change of y = (k x)^p over a step dx from the x that gives the current y.
% y: current value
% dx: step in x
% k: rate
% p: exponent
% dy: change in y

if dx == 0
    dy = 0;
elseif y <= 0
    dy = (k*dx)^p;
else
    x = y^(1/p)/k; % x: x that gives y

    % Slope of (k x)^p:
    dy = p*(k*x)^p/x*dx;
end

end


function dy = dSigmoid(y, dx, yInf, k, p)
% dSigmoid: change of y = 2 yInf (1/2 - 1/(1 + exp((k x)^p))) over a step dx from the x that gives the current y.
% y: current value
% dx: step in x
% yInf: limit of y
% k: rate
% p: exponent
% dy: change in y, zero at or beyond the limit

if dx == 0
    dy = 0;
elseif y <= 0
    dy = 2*yInf*(0.5 - 1/(1 + exp((k*dx)^p)));
else
    e = -(2*yInf/(y - yInf)) - 1; % e: exp((k x)^p) at the x that gives y
    if e <= 1
        dy = 0;
    else
        x = (1/k)*log(e)^(1/p); % x: x that gives y
        z = (k*x)^p; % z: (k x)^p

        % Slope of the sigmoid at x:
        dy = 2*yInf*p*exp(z)*z/(x*(exp(z) + 1)^2)*dx;
    end
end

end
