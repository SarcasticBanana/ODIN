function ps = plotScale(out, cfg)
% plotScale: axis units and averaging choices that keep the demo figures readable at any run length.
% out: state history struct from propagate.m
% cfg: parameter struct from config.m
% ps: struct of plotting settings consumed by main_demo.m

n = numel(out.t); % n: number of output samples
ps.nOrb = out.t(end)/out.T; % ps.nOrb: length of the run [orbits]

% Time unit chosen so the axis reads in tens rather than thousands:
if out.t(end) < 4*3600
    ps.tDiv = 60; ps.tName = 'min';
elseif out.t(end) < 4*86400
    ps.tDiv = 3600; ps.tName = 'h';
elseif out.t(end) < 400*86400
    ps.tDiv = 86400; ps.tName = 'days';
else
    ps.tDiv = 365.25*86400; ps.tName = 'years';
end

ps.t = out.t / ps.tDiv; % ps.t: 1 x N output times in the chosen unit

% An instantaneous trace stays legible while each orbit still spans a few pixels:
ps.dense = ps.nOrb > 30; % ps.dense: true when the full span has to be shown as an average

ps.wOrb = min(max(3, round(out.T/cfg.ode.dtOut)), n); % ps.wOrb: one orbital period [samples]

% Averaging window, at least five orbits so the output grid ripple averages out
% and at most a hundredth of the run so the seasonal trend survives:
ps.nAvg = max(5, round(ps.nOrb/100)); % ps.nAvg: width of the averaging window [orbits]
ps.w = min(max(3, round(ps.nAvg*out.T/cfg.ode.dtOut)), n); % ps.w: averaging window [samples]

% Averaged traces are trimmed by half a window at each end, where the window is short:
if ps.dense
    edge = min(round(ps.w/2), floor((n-1)/2)); % edge: samples discarded at each end
    ps.full = (1 + edge) : (n - edge);
else
    ps.full = 1 : n;
end
% ps.full: samples an averaged trace may be drawn over

% Smooth traces are subsampled for drawing, which instantaneous traces must not be:
ps.dec = unique(round(linspace(ps.full(1), ps.full(end), min(numel(ps.full), 5000))));
% ps.dec: indices of the samples actually drawn for smooth traces

ps.zoom = out.t <= min(3*out.T, out.t(end)); % ps.zoom: mask of the leading orbits drawn instantaneously
ps.zt = out.t(ps.zoom)/60; % ps.zt: zoom panel times [min]

end
