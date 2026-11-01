# ================================================================================
MODULAR ANGULAR RATE PIPELINE

DESCRIPTION:
This software pipeline processes asynchronous event-camera telemetry to detect, track, and characterize faint, unresolved orbital debris (e.g., 1 to 3 cm) moving across the sky. It extracts the kinematic angular rate (arcsec/sec) of the target by utilizing an optimal matched filter and a two-pass dynamic time-binning engine to eliminate motion smear.

The architecture is fully modularized. It separates the time-stationary noise cleanup and data ingestion from the kinematic solver, allowing the foundational event-camera processing functions to be seamlessly reused for downstream size class and spin rate estimation models. This pipeline directly supports the payload team's point-source detection modeling and the mission requirement to analyze the 600 to 1,000 km altitude shell.

# ================================================================================
I. REQUIRED FILES & DEPENDENCIES

To execute this pipeline, ensure the following modular files are present in your MATLAB path:

1. AngularRate_Code.m         : The primary execution script and pipeline coordinator.
2. load_event_data.m          : Ingests, formats, and deduplicates the raw hardware telemetry.
3. build_steady_source_mask.m : Identifies and removes stationary noise (stars and hot pixels).
4. init_matched_filter.m      : Initializes the optical stencil and packages tracking configurations.
5. mf_psf.m                   : Generates the 2D optical Point Spread Function (PSF).
7. track_bins.m               : The core tracking engine that groups events and applies the matched filter.
8. robust_line_fit.m          : The kinematic solver utilizing a Theil-Sen estimator and sigma clipping.

# ================================================================================
II. REQUIRED USER INPUTS (Configured in AngularRate_Code_6.m)

-- File Data --

* fname         : Filepath to the raw event data. Expects a text file where columns represent x, y, polarity, and timestamp.
* zeroIndexed   : Boolean. Set to true if the event camera hardware outputs pixel coordinates starting at 0. MATLAB requires 1-based indexing, so this adds +1 to all coordinates if true.

-- Sensor & Optics Settings --

* nx, ny        : Sensor array resolution in pixels (e.g., 1280 and 720).
* pixPitch_m    : Physical width of a single pixel in meters (e.g., 4.86e-6).
* focalLen_m    : Telescope focal length in meters (e.g., 0.6).
* nw            : Matched-filter template size in pixels (e.g., 21 creates a 21x21 grid).

-- Steady-Source (Noise) Rejection --

* nSeg          : Number of equal time chunks to divide the recording into for noise profiling.
* segFrac       : Fractional threshold (0.0 to 1.0). If a pixel fires in this fraction of segments (e.g., 0.6 = 60%), it is flagged as steady-state noise.
* maskPad       : Dilation padding. Expands the deletion zone by this many pixels around a noisy pixel to prevent bright edges from leaking into the track.

-- Detection Engine Settings --

* faBudget      : Allowed False Alarm Budget (e.g., 0.01 = 1%). Drives the Poisson threshold calculation for spatiotemporal autodetection.
* boxHalf       : Size of the local raw-event counting box used for autodetection around the peak (e.g., 2 creates a 5x5 pixel bounding box).
* cenHalf       : Size of the sub-pixel centroiding window (e.g., 3 creates a 7x7 pixel window).

-- Angular Rate Kinematic Settings --

* dtCoarse      : Shutter speed (seconds) for the Pass 1 coarse tracker (e.g., 0.05).
* pxPerBin      : Target motion per time bin (pixels). Used in Pass 2 to dynamically calculate the perfect integration time to prevent motion smear (e.g., 3).
* clipSigma     : Sigma-clipping threshold for the robust line fit (e.g., 3 standard deviations).
* minResid      : Minimum residual distance floor for the line fit in pixels (e.g., 1.5).

# ================================================================================
III. PIPELINE ARCHITECTURE (Function Breakdown)

* load_event_data.m
Ingests the raw telemetry and strips out identical duplicate events caused by hardware glitches. It rebases the microsecond timestamps so the first recorded event establishes t = 0 seconds.
* build_steady_source_mask.m
Performs time-stationary cleanup by dividing the recording into chunks. Defective hot pixels and stationary stars fire constantly. This function finds pixels that fire constantly across segments, flags them, dilates the boundaries, and masks them out so they cannot dominate the noise profile.
* mf_psf.m
Mathematically models how the target's light spreads across the sensor's pixel grid based on hardware-specific telescope parameters. It creates the 2D optical Point Spread Function (PSF) used as a matched filter template for centroiding.


* init_matched_filter.m
Calls the PSF generator and uses Singular Value Decomposition (SVD) to split the 2D template into two 1D vectors (gx and gy). Sliding 1D vectors is exponentially faster computationally than running full 2D convolutions.
* track_bins.m
The core event-camera tracking engine. It bins events into discrete time chunks and applies the optical matched filter to amplify real targets over noise. It executes the event autodetection logic where it flags a transient detection if the local event count N >= N_min. It ultimately computes an intensity-weighted sub-pixel centroid for the target.


* robust_line_fit.m
Receives the sub-pixel coordinates and calculates velocity. Standard least squares regression is easily skewed by outliers. This function uses a Theil-Sen estimator (calculating the median of all pairwise slopes) to establish a highly robust baseline, followed by an iterative, sigma-clipped least-squares fit to prune stray detections.
* AngularRate_Code_6.m (Main)
Orchestrates the pipeline. It runs the tracker in two passes: a coarse pass to estimate speed, and a fine pass that dynamically scales the integration time to eliminate motion smear. It calculates the physical angular rate using the rate-of-change approximation omega = IFOV * (c(t2) - c(t1)) / (t2 - t1).



# ================================================================================
IV. PHYSICAL & MATHEMATICAL ASSUMPTIONS

* Unresolved Point Sources: The geometric footprint of the orbital debris is treated as an unresolved point source.


* Point Spread Function (PSF) Dispersion: Because the object is unresolved, its light diffuses into an Airy Disk containing ~84% of the source photons in its central lobe. The payload model assumes only a fraction of the PSF (e.g., 60%) falls directly on the peak pixel. The matched filter leverages this exact known physical shape to calculate sub-pixel centroids.


* Poisson Background Statistics: The autodetection step assumes that uncorrelated background noise events arrive randomly. The threshold N_min ensures that any detected cluster represents a statistically significant anomaly against this background, adhering strictly to the false alarm budget.
* Constant Velocity Kinematics: The line-fitting algorithm strictly assumes a constant-velocity moving-point-source model, defined mathematically as theta(t) = theta_0 + omega * t.


* Event Camera Feasibility: Based on the payload team's sensor feasibility model, the minimum photon rate required to trigger the event camera is estimated at 6,000 photons/sec.


* Foundations of Optical Tracking: The tracking theory embedded in this pipeline aligns with the foundational principles of optical multi-object tracking.