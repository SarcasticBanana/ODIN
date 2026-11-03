function [events] = eventDataStreamFun(x0_init, y0_init, vx, vy, Phi_0, T_sim)
    tic;
    %% 0. RANDOM SEED
    % Fixing the seed makes every run identical. Use 'shuffle' for a new random
    % realization each time.
    seed = "shuffle";
    rng(seed);
    
    %% 1. RUN SETTINGS
    bgMode = 'full noise';          % 'clean'  or  'full noise'   (see header)
    
    nx = 2560;  ny = 2560;      % sensor size in pixels. MUST match nx/ny in your
                               % angular-rate pipeline.
    if nargin<6
        T_sim = min(abs((nx.*(vx>0)-x0_init)./vx),abs((ny.*(vy>0)-y0_init)./vy));              % total simulated time (seconds)
    end
    disp(['Simulating T = ',num2str(T_sim), 's of data']);
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
    lambda_r_full = 10.65;       % Reference photon flux (photons/s) in full-noise
                               % mode. The closer it is to lam_bg, the QUIETER the
                               % background. lambda_r_full = lam_bg gives zero
                               % noise. With a=10, Delta_lam=0.10, lam_bg=5:
                               %   4.5 -> ~12 ev/s/px   4.7 -> ~4.5 ev/s/px
                               %   4.8 -> ~1.1 ev/s/px  4.9 -> ~0.1 ev/s/px
    maxEventsAllowed = 10e15;    % Safety cap: full-noise mode stops with an error
                               % if its event estimate exceeds this.
    
    % --- used by BOTH modes ---
    lam_bg = 12;                % Sky background light on each pixel (photons/s)
    
    %% 3. EVENT CIRCUIT PARAMETERS  (the "leaky bucket")
    a         = 2.2;          % Leak (voltage decay) rate, 1/s. Bigger = the
                               % bucket empties faster = shorter tail and (in full
                               % noise mode) more background events.
    b         = 100.0;         % Gain (how much each photon splashes in). It
                               % scales the thresholds equally, so it does not
                               % change the event rate for steady light.
    Delta_lam = 0.13;          % Relative brightness change that triggers an
                               % event (0.10 = 10%). Bigger = fewer events.
    %% 4. TARGET & OPTICS - Validation
    % x0_init = 0;  y0_init = 0;       % Starting centroid of the target (pixels)
    % vx = 600;  vy = 600; % Target velocity (pixels/second). Negative vx =
                               % moving left. Image y points DOWN, so positive vy
                               % moves the target down the screen.
    t_delay = 0.0;
    x0 = x0_init - (vx * t_delay);
    y0 = y0_init - (vy * t_delay);
    
    % Phi_0     = 12000;         % Target's average total photon flux (photons/s),
                               % summed over all the pixels its blur covers
    mod_depth = 0.0;           % Brightness flicker size as a fraction of Phi_0
                               % (0 = steady, 0.5 = +/-50%, 1 = fades to zero)
    f_spin    = 10;            % Flicker frequency in Hz (ignored if mod_depth=0)
    sigma_psf = 1.15;           % Blur (Gaussian PSF standard deviation) in pixels
            
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
    if nargin<6
        writematrix(events_out, sprintf("model/Payload Models/Event Data Simulation Model/census/synthetic_debris_events_x%i_y%i_vx%.2f_vy%.2f_tsim%.2f.txt",x0,y0_init,vx,vy,T_sim), 'Delimiter', 'space');
    else
        writematrix(events_out, sprintf("model/Payload Models/Event Data Simulation Model/censusNoise/synthetic_debris_events_x%i_y%i_vx%.2f_vy%.2f_tsim%.2f.txt",x0,y0_init,vx,vy,T_sim), 'Delimiter', 'space');
    end
    % fprintf('Wrote synthetic_debris_events.txt  (columns: x y pol t_us)\n');
end

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