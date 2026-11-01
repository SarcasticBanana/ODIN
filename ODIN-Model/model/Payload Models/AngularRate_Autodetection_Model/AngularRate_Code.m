%% MODULAR_ANGULAR_RATE_PIPELINE
% Extracts the kinematic angular rate of a target.
% Built on a reusable foundational framework for event camera processing.
% REQUIRES: track_bins.m and matchedfilter_psf.m on the MATLAB path.

clear; clc; close all;
tStart = tic;

%% 1. USER INPUTS
 %fname = "C:\Users\HuntR\Aero Undergrad\FA26\Senior Projects\ONEWEB_0023_2024-11-15_18-50-16_cut.txt";
% % --- Sensor & Optics Settings ---
 nx = 2560;  ny = 2560; % pixels 
 zeroIndexed = false; 

fname = "C:\Users\HuntR\Aero Undergrad\FA26\Senior Projects\synthetic_debris_events.txt";
%nx = 1280;  ny = 720;        % Match the grid size from your simulator (or 1280x720
%zeroIndexed = false;      % Simulator outputs 1-based MATLAB coordinates (1 to nx)

if isstring(fname); fname = convertStringsToChars(fname); end
pixPitch_m = 4.86e-6; % meters            
focalLen_m = 0.475;     % meters             
IFOV_rad   = 2*atan(pixPitch_m/(2*focalLen_m));            
nw = 21;                           

% --- Steady-Source (Noise) Rejection ---
nSeg    = 10;                      
segFrac = 0.6;                     
maskPad = 2;                       

% --- Detection Settings ---
faBudget = 0.01;                   
boxHalf  = 2;                      
cenHalf  = 3;                      

% --- Angular Rate Specific Settings ---
dtCoarse = 0.05;                   
pxPerBin = 3;                      
clipSigma = 3;                     
minResid  = 1.5;                   

%% 2. EXECUTE REUSABLE FOUNDATION
% 1. Load and format raw asynchronous events
[t, x, y, T_duration] = load_event_data(fname, nx, ny, zeroIndexed);

% 2. Build the spatial noise mask (removes stationary stars and hot pixels)
steadyMask = build_steady_source_mask(t, x, y, T_duration, nx, ny, nSeg, segFrac, maskPad);

% 3. Generate the optical matched-filter stencil and package it into a config struct
cfg = init_matched_filter(nw, nx, ny, steadyMask, boxHalf, cenHalf, faBudget);

%% 3. EXECUTE UNIFIED ANGULAR RATE PIPELINE
% Runs the coarse tracker, fine tracker, kinematic solver, console report, and plot generation
results = analyze_angular_rate(t, x, y, nx, ny, steadyMask, cfg, IFOV_rad, dtCoarse, pxPerBin, clipSigma, minResid, tStart);


%% =======================================================================
%               REUSABLE FOUNDATION FUNCTIONS
% ========================================================================
%% =======================================================================
%               REUSABLE FOUNDATION FUNCTIONS
% ========================================================================
% 
% function [t, x, y, T_duration] = load_event_data(fname, nx, ny, zeroIndexed)
% 
% function steadyMask = build_steady_source_mask(t, x, y, T, nx, ny, nSeg, segFrac, maskPad)
% 
% function cfg = init_matched_filter(nw, nx, ny, steadyMask, boxHalf, cenHalf, faBudget)

%% =======================================================================
%               UNIFIED ANGULAR RATE FUNCTION
% ========================================================================

function results = analyze_angular_rate(t, x, y, nx, ny, steadyMask, cfg, IFOV_rad, dtCoarse, pxPerBin, clipSigma, minResid, tStart)
    
    % ==========================================
    % PHASE 1: KINEMATIC SOLVING
    % ==========================================
    
    % Pass 1: Run the tracker at a coarse, fixed shutter speed (dtCoarse) 
    trk1 = track_bins(t, x, y, dtCoarse, cfg);
    assert(size(trk1,1) >= 5, 'Pass 1 found fewer than 5 detections.');
    [px1, py1, ~] = robust_line_fit(trk1(:,1), trk1(:,2), trk1(:,3), clipSigma, minResid);
    v1 = hypot(px1(1), py1(1));
    
    % Pass 2: Dynamically recalculate the shutter speed (dtFine)
    dtFine = min(dtCoarse, max(0.001, round(pxPerBin/v1, 3)));
    trk = track_bins(t, x, y, dtFine, cfg);
    assert(size(trk,1) >= 5, 'Pass 2 found fewer than 5 detections.');
    
    % Extract the final, highly-accurate timeline and sub-pixel coordinates
    tt = trk(:,1);  xs = trk(:,2);  ys = trk(:,3);
    [px, py, keep, res] = robust_line_fit(tt, xs, ys, clipSigma, minResid);
    
    % Calculate the standard error of the fitted line slope
    tk  = tt(keep);  nk = numel(tk);
    Stt = sum((tk-mean(tk)).^2);
    sx  = sqrt(sum((xs(keep)-polyval(px,tk)).^2)/(nk-2))/sqrt(Stt);   
    sy  = sqrt(sum((ys(keep)-polyval(py,tk)).^2)/(nk-2))/sqrt(Stt);
    
    % Convert the final pixel velocities into real-world Angular Rates
    vx = px(1);  vy = py(1);  vmag = hypot(vx,vy);                    
    rad2as = 180/pi*3600;
    
    % Package deliverables into a single results structure
    results.w_x = IFOV_rad*vx*rad2as;   
    results.w_y = IFOV_rad*vy*rad2as;             
    results.w   = IFOV_rad*vmag*rad2as;
    results.sw  = IFOV_rad*sqrt((vx*sx)^2 + (vy*sy)^2)/vmag*rad2as;
    results.heading = atan2d(vy,vx);           
    results.vx = vx; results.vy = vy; results.vmag = vmag;
    results.sx = sx; results.sy = sy;
    results.dtFine = dtFine; results.IFOV = IFOV_rad*rad2as;
    
    ik = find(keep);  i0 = ik(1);  i1 = ik(end);
    results.w2pt = IFOV_rad*hypot(xs(i1)-xs(i0), ys(i1)-ys(i0))/(tt(i1)-tt(i0))*rad2as;
    
    % ==========================================
    % PHASE 2: CONSOLE REPORTING
    % ==========================================
    fprintf('\n=========== ANGULAR RATE RESULT ===========\n');
    fprintf('Detections: %d   Inliers: %d (%.1f%%)   rms residual: %.2f px\n', ...
        size(trk,1), nk, 100*nk/size(trk,1), sqrt(mean(res(keep).^2)));
    fprintf('Transit: %.2f to %.2f s (%.2f s)   bin dt = %.3f s   motion = %.1f px/bin\n', ...
        min(tk), max(tk), max(tk)-min(tk), results.dtFine, results.vmag*results.dtFine);
    fprintf('Pixel velocity : vx = %+.2f +/- %.2f px/s,  vy = %+.2f +/- %.2f px/s\n', results.vx, results.sx, results.vy, results.sy);
    fprintf('Angular rate   : wx = %+.1f arcsec/s,  wy = %+.1f arcsec/s\n', results.w_x, results.w_y);
    fprintf('|omega| = %.1f +/- %.1f arcsec/s = %.4f deg/s (statistical, 1-sigma)\n', ...
        results.w, results.sw, results.w/3600);
    fprintf('Heading = %.1f deg (from +x toward +y; image y points down)\n', results.heading);
    fprintf('Two-point slide formula (first/last inlier): %.1f arcsec/s\n', results.w2pt);
    
    if results.vmag*results.dtFine > 5
        warning('Motion per bin is %.1f px (> 5): expect smear. Lower pxPerBin or dtCoarse.', results.vmag*results.dtFine);
    end
    fprintf('Elapsed: %.1f s\n', toc(tStart));

    % ==========================================
    % PHASE 3: PLOT GENERATION
    % ==========================================
    Nimg = accumarray([y x],1,[ny nx]);
    Nimg(steadyMask) = 0;              
    tl = [min(tk) max(tk)];

    figure('Position',[50 50 1300 800]);
    
    subplot(2,2,1);
    imagesc(Nimg); axis image; colormap(gca,hot(256)); caxis([0 3]); hold on
    plot(xs(keep),ys(keep),'g.','MarkerSize',5);
    plot(xs(~keep),ys(~keep),'rx');
    plot(polyval(px,tl),polyval(py,tl),'c-','LineWidth',1);
    xlabel('x [px]'); ylabel('y [px]');
    title(sprintf('Transient events and track:  |\\omega| = %.1f arcsec/s', results.w));
    legend('inlier','rejected','fit','Location','southoutside','Orientation','horizontal');

    subplot(2,2,2);
    plot(tt(keep),xs(keep),'g.'); hold on; plot(tt(~keep),xs(~keep),'rx');
    plot(tl,polyval(px,tl),'k-'); grid on
    xlabel('Time [s]'); ylabel('x [px]'); title(sprintf('v_x = %+.1f px/s', results.vx));

    subplot(2,2,3);
    plot(tt(keep),ys(keep),'g.'); hold on; plot(tt(~keep),ys(~keep),'rx');
    plot(tl,polyval(py,tl),'k-'); grid on
    xlabel('Time [s]'); ylabel('y [px]'); title(sprintf('v_y = %+.1f px/s', results.vy));

    subplot(2,2,4);
    plot(tt(keep),res(keep),'g.'); hold on; plot(tt(~keep),res(~keep),'rx'); grid on
    ylim([0 max(3,1.5*max(res(keep)))]);
    xlabel('Time [s]'); ylabel('Distance from fit [px]'); title('Residuals');
end
