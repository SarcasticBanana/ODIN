%% EVENT_DATA_STREAM_SIMULATION_V3
% Simulates the output of an event camera watching one bright point source
% move across the sky. Based on Benson & Holzinger (2022), Eqs. 1-15.
%
% OUTPUT FILE:  synthetic_debris_events.txt
%               one event per line:   x   y   polarity   time_in_microseconds
%
% -------------------------------------------------------------------------
% THE BIG IDEA (read this first)
% -------------------------------------------------------------------------
% An event camera does not take pictures. Every pixel is its own little
% "leaky bucket":
%   * Each photon that hits the pixel pours a splash of water into the bucket
%     (the voltage V jumps up).
%   * The bucket leaks constantly, so V drifts back toward 0 (decay rate 'a').
%   * The pixel remembers a "last checked" level, V_ref.
%   * If V rises  c_pos above V_ref  -> the pixel fires a  +1 event.
%     If V falls  |c_neg| below V_ref -> the pixel fires a  -1 event.
%     After any event, V_ref is reset to the current V.
% A target passing over a pixel makes V shoot up (+ events: the leading
% "bow shock"), then V leaks back down afterward (- events: the "tail").
%
% -------------------------------------------------------------------------
% TWO BACKGROUND (NOISE) MODES
% -------------------------------------------------------------------------
%   'clean'      Target-only. The reference flux lambda_r is set equal to the
%                sky flux lam_bg, so sky photons add ZERO voltage and cause
%                ZERO events. Only pixels near the track are simulated (fast).
%                Optional sparse random noise (noiseRate) is sprinkled on top.
%
%   'full noise' The sky itself makes every pixel fire, because lambda_r is
%                set below lam_bg (see lambda_r_full). Every pixel on the
%                sensor must be simulated (slow, many events). The amount of
%                noise is tuned with lambda_r_full. noiseRate is NOT used.

% ANALYSIS / CALIBRATION PLOTS: run Event_Data_Diagnostics.m on the output file.

clear; clc; close all; tic;

%% 0. RANDOM SEED
% Fixing the seed makes every run identical. Use 'shuffle' for a new random
% realization each time.
seed = 1;
rng(seed);

%% 1. RUN SETTINGS
bgMode = 'clean';          % 'clean'  or  'full noise'   (see header)

nx = 1280;  ny = 720;      % sensor size in pixels. MUST match nx/ny in your
                           % angular-rate pipeline.
T_sim = 2.09;              % total simulated time (seconds)
dt    = 1e-3;              % time step used to LOOK UP the target's light on
                           % each pixel (seconds). Photon arrival times are
                           % NOT limited to this grid. Keep the target's
                           % motion per step (v*dt) below ~0.3 px.

polarity01 = true;         % true : write polarity as 0/1 (Prophesee style)
                           % false: write polarity as -1/+1

%% 2. NOISE INPUTS  (everything that controls how noisy the data is)
% --- used by 'clean' mode only ---
noiseRate = 0.10;          % Random "salt" noise added on top of the target
                           % events, in events / second / pixel. Positions,
                           % times and polarities are uniformly random.
                           % Calibrate this to your real data's median rate.

% --- used by 'full noise' mode only ---
lambda_r_full = 4.8;       % Reference photon flux (photons/s) in full-noise
                           % mode. The closer it is to lam_bg, the QUIETER the
                           % background. lambda_r_full = lam_bg gives zero
                           % noise. With a=10, Delta_lam=0.10, lam_bg=5:
                           %   4.5 -> ~12 ev/s/px   4.7 -> ~4.5 ev/s/px
                           %   4.8 -> ~1.1 ev/s/px  4.9 -> ~0.1 ev/s/px
maxEventsAllowed = 3e7;    % Safety cap: full-noise mode stops with an error
                           % if its event estimate exceeds this.

% --- used by BOTH modes ---
lam_bg = 5;                % Sky background light on each pixel (photons/s)

%% 3. EVENT CIRCUIT PARAMETERS  (the "leaky bucket")
a         = 10.0;          % Leak (voltage decay) rate, 1/s. Bigger = the
                           % bucket empties faster = shorter tail and (in full
                           % noise mode) more background events.
b         = 100.0;         % Gain (how much each photon splashes in). It
                           % scales the thresholds equally, so it does not
                           % change the event rate for steady light.
Delta_lam = 0.10;          % Relative brightness change that triggers an
                           % event (0.10 = 10%). Bigger = fewer events.

%% 4. TARGET & OPTICS
x0 = 1280;  y0 = 90;       % Starting centroid of the target (pixels)
vx = -612.44;  vy = 43.06; % Target velocity (pixels/second). Negative vx =
                           % moving left. Image y points DOWN, so positive vy
                           % moves the target down the screen.
Phi_0     = 12000;         % Target's average total photon flux (photons/s),
                           % summed over all the pixels its blur covers
mod_depth = 0.0;           % Brightness flicker size as a fraction of Phi_0
                           % (0 = steady, 0.5 = +/-50%, 1 = fades to zero)
f_spin    = 10;            % Flicker frequency in Hz (ignored if mod_depth=0)
sigma_psf = 1.5;           % Blur (Gaussian PSF standard deviation) in pixels

% Optics, used ONLY to print the ground-truth angular rate at the end
pixPitch_m = 4.86e-6;      % pixel size (meters)
focalLen_m = 0.6;          % telescope focal length (meters)

%% 5. DERIVED QUANTITIES
fullNoise = strcmpi(bgMode, 'full noise');
if ~fullNoise && ~strcmpi(bgMode, 'clean')
    error('bgMode must be ''clean'' or ''full noise''.');
end

% Reference flux: sets where V = 0 sits. Clean mode puts it exactly at the
% sky level so the sky is invisible; full-noise mode puts it lower.
if fullNoise, lambda_r = lambda_r_full; else, lambda_r = lam_bg; end

% Event thresholds (Eq. 5). A relative brightness change of Delta_lam at
% steady state moves the voltage by exactly this much. c_neg is negative.
c_pos = (b/a)*log(1 + Delta_lam);
c_neg = (b/a)*log(1 - Delta_lam);

% Bundle circuit constants to hand to the pixel function
p = struct('a',a, 'b',b, 'lambda_r',lambda_r, 'c_pos',c_pos, 'c_neg',c_neg);

%% 6. TARGET MOTION AND LIGHT OVER TIME
t_vec  = (0:dt:T_sim)';                       % time grid (column vector)
x_bar  = x0 + vx*t_vec;                       % target centroid x at each time
y_bar  = y0 + vy*t_vec;                       % target centroid y at each time
Phi_t  = Phi_0*(1 + mod_depth*sin(2*pi*f_spin*t_vec));  % target brightness vs time
inv_s2 = 1/(sqrt(2)*sigma_psf);               % constant used inside erf()

% Warn if the blur would be clipped by the sensor edge
m = ceil(5*sigma_psf) + 1;                    % margin: pixels beyond ~5 sigma see no target
if min(x_bar) < 1+m || max(x_bar) > nx-m || min(y_bar) < 1+m || max(y_bar) > ny-m
    warning('Track leaves the sensor or passes within %d px of an edge.', m);
end

%% 7. CHOOSE WHICH PIXELS TO SIMULATE
if fullNoise
    % Every pixel fires from sky photons, so every pixel must be simulated
    xList = 1:nx;  yList = 1:ny;
else
    % Clean mode: pixels away from the track only see steady sky light, which
    % makes no events, so skip them. Simulate just the box around the track.
    xList = max(1,floor(min(x_bar))-m) : min(nx,ceil(max(x_bar))+m);
    yList = max(1,floor(min(y_bar))-m) : min(ny,ceil(max(y_bar))+m);
end

% Full-noise safety check: rough estimate of how many events will be made.
% Each sky photon causes 1 "+" event, then the voltage decays back down in
% steps of |c_neg|, causing about h/|c_neg| "-" events. (A conservative
% over-estimate when lambda_r is close to lam_bg.)
if fullNoise
    h_bg        = b*(log(lam_bg) - log(lambda_r))/lam_bg;   % voltage jump per sky photon
    evPerPhoton = 1 + max(0, h_bg)/abs(c_neg);
    estEvents   = numel(xList)*numel(yList)*T_sim*lam_bg*evPerPhoton;
    fprintf('Estimated background events: %.2e\n', estEvents);
    if estEvents > maxEventsAllowed
        error(['Full-noise mode would create ~%.1e events. Reduce nx/ny, T_sim, ' ...
               'or a, raise lambda_r_full, or raise maxEventsAllowed.'], estEvents);
    end
end

fprintf('%s mode: simulating %d x %d = %d pixels of a %dx%d sensor\n', ...
    bgMode, numel(xList), numel(yList), numel(xList)*numel(yList), nx, ny);

%% 8. SIMULATE EVERY SELECTED PIXEL
% The target is a 2-D Gaussian blur. The fraction of its light landing in a
% pixel is a product of an x part and a y part, each computed with erf()
% (Eq. 15). Fy depends only on the pixel row, so compute it once for all rows.
%   Fy size: (number of time steps) x (number of rows)
Fy = 0.5*(erf((yList + 0.5 - y_bar)*inv_s2) - erf((yList - 0.5 - y_bar)*inv_s2));

events    = zeros(1e6, 4);     % columns: x, y, polarity(-1/+1), time(s)
evCount   = 0;                 % how many rows of 'events' are filled
nCols     = numel(xList);
nextPrint = 0.1;               % next progress-print threshold (10%, 20%, ...)

for ci = 1:nCols
    ix = xList(ci);
    % Fraction of the target's light in this pixel column, at every time step
    Fx = 0.5*(erf((ix + 0.5 - x_bar)*inv_s2) - erf((ix - 0.5 - x_bar)*inv_s2));

    for cj = 1:numel(yList)
        iy = yList(cj);

        % Photon rate (photons/s) hitting THIS pixel at every time step:
        % steady sky + the target's moving, blurred light
        lam_t = lam_bg + Phi_t .* Fx .* Fy(:,cj);

        % Run the leaky-bucket simulation for this one pixel
        [tEv, pEv] = pixel_stochastic(lam_t, dt, T_sim, p);

        % Append this pixel's events to the big list, doubling its size
        % whenever it fills up (much faster than growing it each time)
        k = numel(tEv);
        if k > 0
            need = evCount + k;
            if need > size(events,1)
                events(max(2*size(events,1), need), 4) = 0;
            end
            events(evCount+1:need,:) = [repmat(ix,k,1), repmat(iy,k,1), pEv, tEv];
            evCount = need;
        end
    end

    if ci/nCols >= nextPrint
        fprintf('  %3.0f%%  (%d events so far)\n', 100*ci/nCols, evCount);
        nextPrint = nextPrint + 0.1;
    end
end
events  = events(1:evCount,:);   % drop the unused pre-allocated rows
nTarget = evCount;               % events produced by the circuit simulation

%% 9. ADD SPARSE RANDOM NOISE (clean mode only)
% Total noise events ~ rate x pixels x time. A Poisson count is approximated
% with a normal distribution (needs no toolbox). Each noise event gets a
% uniformly random pixel, polarity and time.
if ~fullNoise && noiseRate > 0
    mu    = noiseRate*nx*ny*T_sim;
    nN    = max(0, round(mu + sqrt(mu)*randn));
    noise = [randi(nx,nN,1), randi(ny,nN,1), 2*(rand(nN,1)>0.5)-1, rand(nN,1)*T_sim];
    events = [events; noise];
    fprintf('Added %d sparse noise events (%.3g events/s/pixel)\n', nN, noiseRate);
end

if isempty(events), warning('No events triggered!'); return; end

%% 10. SORT BY TIME AND SAVE
[~, sortIdx] = sort(events(:,4));
events  = events(sortIdx,:);
pol_out = events(:,3);
if polarity01, pol_out = (pol_out + 1)/2; end       % -1/+1  ->  0/1
events_out = [events(:,1), events(:,2), pol_out, round(events(:,4)*1e6)];  % time in microseconds

fprintf('Finished in %.1f s. Total events: %d (circuit %d; +: %d, -: %d)\n', ...
    toc, size(events_out,1), nTarget, sum(events(:,3)>0), sum(events(:,3)<0));
writematrix(events_out, 'synthetic_debris_events.txt', 'Delimiter', 'space');
fprintf('Wrote synthetic_debris_events.txt  (columns: x y pol t_us)\n');

%% 11. GROUND TRUTH  (the answer your angular-rate pipeline should reproduce)
IFOV_as = 2*atan(pixPitch_m/(2*focalLen_m))*180/pi*3600;   % arcsec per pixel
v_true  = hypot(vx, vy);                                   % pixel speed (px/s)
fprintf('\n---- GROUND TRUTH ----\n');
fprintf('vx = %+.2f px/s, vy = %+.2f px/s, |v| = %.2f px/s\n', vx, vy, v_true);
fprintf('IFOV = %.3f arcsec/px  ->  |omega| = %.2f arcsec/s (%.4f deg/s)\n', ...
    IFOV_as, IFOV_as*v_true, IFOV_as*v_true/3600);
fprintf('heading = %.2f deg, transit %.2f s, %.1f px\n', atan2d(vy,vx), T_sim, v_true*T_sim);
fprintf('\nFor Event_Data_Diagnostics.m, set:\n');
fprintf('  trackX0 = %g;  trackY0 = %g;  trackVx = %g;  trackVy = %g;  trackT0 = 0;\n', x0, y0, vx, vy);

%% 12. PLOTS  (two quick-look figures)
% The calibration plots (windowed snapshot, event rate vs time, offsets from
% the track) live in Event_Data_Diagnostics.m so the SAME code can be run on
% real and synthetic data. Here are just two quick checks that the target
% moved the way you asked.
nTb       = 200;           % number of time bins for the space-time plots
showTruth = true;          % overlay the true target path as a dashed line

tb    = min(floor(events(:,4)/T_sim*nTb) + 1, nTb);   % time-bin index of each event
tAxis = ((1:nTb)-0.5)*T_sim/nTb;                      % center time of each bin (s)

% ---- FIGURE 1: x vs time and y vs time -----------------------------------
% Each plot counts events per (time bin, position). The target is a straight
% line; its slope is the pixel velocity. Everything else is background.
cntX = accumarray([tb, events(:,1)], 1, [nTb nx]);
cntY = accumarray([tb, events(:,2)], 1, [nTb ny]);

figure('Color','w','Name','Space-time density','Position',[100 100 1150 450]);
subplot(1,2,1);
imagesc(tAxis, 1:nx, log10(1+cntX')); axis xy; colorbar; hold on;
if showTruth, plot(t_vec, x_bar, 'w--'); end
ylim([1 nx]);
xlabel('Time (s)'); ylabel('x [px]'); title('x vs time, log_{10}(1+N)');

subplot(1,2,2);
imagesc(tAxis, 1:ny, log10(1+cntY')); axis xy; colorbar; hold on;
if showTruth, plot(t_vec, y_bar, 'w--'); end
ylim([max(1, min(y_bar)-40), min(ny, max(y_bar)+40)]);   % zoom on the rows the target crosses
xlabel('Time (s)'); ylabel('y [px]'); title('y vs time, log_{10}(1+N)');

% ---- FIGURE 2: whole-run accumulated event map ---------------------------
% Every event over the entire run, drawn on the sensor. The target leaves a
% bright streak; the faint speckle is background noise.
Nimg = accumarray([events(:,2), events(:,1)], 1, [ny nx]);

figure('Color','w','Name','Accumulated events','Position',[100 100 1000 600]);
imagesc(log10(1+Nimg)); axis image; colormap(gca, hot(256)); colorbar; hold on;
if showTruth, plot(x_bar, y_bar, 'c--'); end
xlim([1 nx]); ylim([1 ny]);
xlabel('x [px]'); ylabel('y [px]');
title('All events accumulated over the run, log_{10}(1+N)  (cyan = true path)');


%% =======================================================================
%                              LOCAL FUNCTION
% ========================================================================
function [tEv, pEv] = pixel_stochastic(lam_t, dt, T_sim, p)
% Simulates ONE pixel photon-by-photon and returns its events.
%   lam_t : photon rate hitting the pixel at each time step (photons/s)
%   dt    : spacing of lam_t samples (s)
%   T_sim : total duration (s)
%   p     : circuit constants (a, b, lambda_r, c_pos, c_neg)
% Returns:
%   tEv   : event times (s)         pEv : event polarities (+1 / -1)

    a = p.a;  b = p.b;  lr = p.lambda_r;  cp = p.c_pos;  cn = p.c_neg;
    nSteps = numel(lam_t);

    % Start the bucket at its steady-state level for the starting light (Eq. 9)
    % so the pixel does not make a burst of fake events at t = 0.
    V     = (b/a)*log(lam_t(1)/lr);
    V_ref = V;                                   % the "last checked" level
    lam_max = max(lam_t);                        % brightest the pixel ever gets

    % ---- STEP 1: generate random photon arrival times -----------------------
    % The light level changes with time, so use "thinning": first create
    % arrivals at the CONSTANT maximum rate, then randomly throw some away
    % so the survivors follow the true time-varying rate (Eq. 8, Sec. 2.2).
    n_est  = ceil(1.2*lam_max*T_sim + 10*sqrt(lam_max*T_sim) + 20);  % how many to draw
    T_cand = cumsum(-log(1 - rand(n_est,1))/lam_max);   % gaps are exponential -> cumulative sum = arrival times
    while T_cand(end) < T_sim                           % drew too few? draw more
        T_cand = [T_cand; T_cand(end) + cumsum(-log(1 - rand(n_est,1))/lam_max)]; %#ok<AGROW>
    end
    T_cand = T_cand(T_cand <= T_sim);                   % keep only arrivals inside the run

    if isempty(T_cand)
        T_n = zeros(0,1);  lam_n = zeros(0,1);
    else
        idx  = min(max(round(T_cand/dt) + 1, 1), nSteps);  % time -> index into lam_t
        lamc = lam_t(idx);                                 % true light level at each candidate
        keep = rand(size(T_cand)) <= lamc/lam_max;         % keep with probability lam/lam_max
        T_n  = T_cand(keep);                               % final photon arrival times
        lam_n = lamc(keep);                                % light level at each photon
    end

    % ---- STEP 2: march through the photons, tracking the voltage ------------
    cap = 256;                                   % size of the event buffer (grows by doubling)
    tEv = zeros(cap,1);  pEv = zeros(cap,1);  nEv = 0;
    nPh    = numel(T_n);
    t_curr = 0;                                  % time up to which V is known

    % One loop pass per photon, plus ONE EXTRA pass for the quiet stretch after
    % the last photon (so its leftover decay events are not missed).
    for n = 1:nPh+1
        isTail = (n > nPh);
        if isTail, t_next = T_sim; else, t_next = T_n(n); end

        % -- 2a. Events caused purely by the voltage LEAKING before the next photon.
        % Between photons V decays smoothly, so we can solve exactly for the
        % moment it reaches the next threshold and, if that moment is before
        % the next photon, record an event and repeat.
        for guard = 1:2000
            if V > 0
                tgt = V_ref + cn;                          % lower threshold level
                if tgt > 0 && V > tgt, pol = -1; else, break; end
            elseif V < 0
                tgt = V_ref + cp;                          % upper threshold (V climbing toward 0)
                if tgt < 0 && V < tgt, pol = +1; else, break; end
            else
                break;
            end
            t_cross = t_curr + log(V/tgt)/a;               % exact time V reaches the threshold
            if t_cross < t_next
                V = tgt;  V_ref = tgt;  t_curr = t_cross;  % event: reset reference
                nEv = nEv + 1;
                if nEv > cap, cap = 2*cap; tEv(cap) = 0; pEv(cap) = 0; end
                tEv(nEv) = t_curr;  pEv(nEv) = pol;
            else
                break;                                     % next threshold not reached before next photon
            end
        end

        % -- 2b. Let V decay the rest of the way to the next photon (or run end)
        V = V*exp(-a*(t_next - t_curr));
        t_curr = t_next;
        if isTail, break; end                              % the extra pass ends here

        % -- 2c. The photon arrives: voltage jumps up (Eq. 3)
        V = V + b*(log(lam_n(n)) - log(lr))/lam_n(n);

        % -- 2d. Did that jump cross a threshold? (at most ONE event per photon)
        d = V - V_ref;
        if d >= cp
            V_ref = V;  pol = +1;
        elseif d <= cn
            V_ref = V;  pol = -1;
        else
            continue;                                      % no event, go to next photon
        end
        nEv = nEv + 1;
        if nEv > cap, cap = 2*cap; tEv(cap) = 0; pEv(cap) = 0; end
        tEv(nEv) = t_curr;  pEv(nEv) = pol;
    end

    tEv = tEv(1:nEv);  pEv = pEv(1:nEv);            % trim unused buffer space
end