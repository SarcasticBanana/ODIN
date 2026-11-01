
function trk = track_bins(t,x,y,dt,cfg)
% TRACK_BINS Generic Event Camera Tracking Engine
% Groups events into time chunks, applies an optical matched filter to amplify 
% real targets, verifies the target using a Poisson false-alarm budget, and 
% calculates an intensity-weighted sub-pixel center of mass.
%
% OUTPUT:
% trk = [t_mid, x_centroid, y_centroid, N_local, N_min]

    bin   = floor(t/dt) + 1;                    % Assign every single event to a specific time bin
    nBins = bin(end);
    cnt   = accumarray(bin,1,[nBins 1]);
    first = [1; 1+cumsum(cnt)];                 

    nFree   = nnz(~cfg.mask);                   % Number of unmasked, usable pixels on the sensor
    boxArea = (2*cfg.boxHalf+1)^2;
    nBox    = nFree/boxArea;                    % How many independent 5x5 testing boxes fit on the screen

    trk = nan(nBins,5);
    for k = 1:nBins
        if cnt(k) < 2, continue; end            % Skip this time bin if fewer than 2 events fired
        
        % Build a 2D image array of all the events that fired during this specific time bin
        ii = first(k):first(k+1)-1;
        N = accumarray([y(ii) x(ii)],1,[cfg.ny cfg.nx]);
        N(cfg.mask) = 0;                        % Zero out the hot pixels and stars
        nEv = sum(N(:));
        if nEv == 0, continue; end

        % Apply the Matched Filter: Slide the 1D X and Y stencils over the 2D event image.
        % This mathematically amplifies clusters that look like a real optical point source.
        MF = conv2(cfg.gy, cfg.gx.', N, 'same');   
        
        % Find the exact pixel coordinate where the stencil matched the data best
        [~,imax] = max(MF(:));
        [iy,ix]  = ind2sub([cfg.ny cfg.nx], imax);

        % Count the raw raw events in a small box around that peak
        ya = max(iy-cfg.boxHalf,1);  yb = min(iy+cfg.boxHalf,cfg.ny);
        xa = max(ix-cfg.boxHalf,1);  xb = min(ix+cfg.boxHalf,cfg.nx);
        nLoc = sum(N(ya:yb,xa:xb),'all');

        % --- Autodetection (Poisson False-Alarm Test) ---
        % Even pure static noise occasionally clusters together. We calculate N_min: the absolute 
        % minimum number of events required in our box to guarantee (with 99% confidence) that 
        % this cluster is real debris and not just a random coincidence of background noise.
        muBg = nEv/nFree*boxArea;               % The average background noise per box
        nMin = 1;
        while nBox*gammainc(muBg,nMin) > cfg.fa && nMin < 1000
            nMin = nMin + 1;
        end
        if nLoc < nMin, continue; end           % If the box event count is lower than N_min, reject it as noise

        % --- Sub-Pixel Centroiding ---
        % We know the peak pixel (ix, iy), but the true center is likely somewhere between pixels.
        % We calculate a weighted average (center of mass) of the event intensities inside a 7x7 box.
        ya = max(iy-cfg.cenHalf,1);  yb = min(iy+cfg.cenHalf,cfg.ny);
        xa = max(ix-cfg.cenHalf,1);  xb = min(ix+cfg.cenHalf,cfg.nx);
        W = N(ya:yb,xa:xb);
        [XX,YY] = meshgrid(xa:xb,ya:yb);
        Itot = sum(W(:));
        % Save [Time_midpoint, X_subpixel, Y_subpixel, Local_events, Required_events]
        trk(k,:) = [(k-0.5)*dt, sum(W(:).*XX(:))/Itot, sum(W(:).*YY(:))/Itot, nLoc, nMin];
    end
    trk = trk(~isnan(trk(:,1)),:); % Strip out empty rows
end

