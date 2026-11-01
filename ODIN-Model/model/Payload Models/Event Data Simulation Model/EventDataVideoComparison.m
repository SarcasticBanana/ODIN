%% REALTIME_RAW_EVENT_STREAM_COMPARISON
% Strictly loads and visually compares two raw event camera data streams 
% side-by-side in true real-time (1.0x wall-clock synchronized).
%
% EXPECTED FILE FORMAT (4 whitespace- or comma-separated numeric columns):
%   x   y   polarity   time_in_microseconds

clear; clc; close all;

%% 1. INPUT FILES & DISPLAY SETTINGS
fileA     = "C:\Users\HuntR\Aero Undergrad\FA26\Senior Projects\ONEWEB_0023_2024-11-15_18-50-16_cut.txt";
titleA    = 'Original recording';
subtitleA = 'Unchanged recording';

%fileB     = "C:\Users\HuntR\Aero Undergrad\FA26\Senior Projects\synthetic_debris_events.txt";
fileB     = "C:\Users\HuntR\VADeR Research\DUBHE_2025-02-26_20-58-47.txt";
titleB    = 'Tuned candidate';
subtitleB = 'Comparison stream';

% Indexing convention for pixel coordinates (x, y):
%   true  -> file uses 0-based indexing (0 to nx-1); converts to 1-based
%   false -> file already uses 1-based indexing (1 to nx)
zeroBasedA = true;
zeroBasedB = false;

% Sensor & Real-Time Replay Parameters
nx             = 1280        % Sensor width (pixels)
ny             = 720;         % Sensor height (pixels)
replayWindow_s = 0.100;       % Time window accumulated per frame (e.g., 0.100 = 100 ms)
playbackSpeed  = 1.0;         % 1.0 = true real-time, 0.5 = half-speed, 2.0 = double-speed
targetFPS      = 60;          % Maximum display refresh rate (Hz)

%% 2. LOAD & PRE-INDEX BOTH EVENT STREAMS
[streamA, durA] = load_and_preindex_stream(fileA, zeroBasedA, nx, ny, titleA);
[streamB, durB] = load_and_preindex_stream(fileB, zeroBasedB, nx, ny, titleB);

if streamA.count == 0 && streamB.count == 0
    error('Neither event file could be loaded. Check file names and paths.');
end

%% 3. LAUNCH REAL-TIME VIDEO COMPARISON
play_realtime_comparison(streamA, durA, titleA, subtitleA, ...
                         streamB, durB, titleB, subtitleB, ...
                         nx, ny, replayWindow_s, playbackSpeed, targetFPS);


%% ========================================================================
%                              LOCAL FUNCTIONS
% =========================================================================

function [s, duration] = load_and_preindex_stream(filePath, isZeroBased, nx, ny, label)
% Loads an event file and pre-computes 1D linear pixel indices so the 
% real-time render loop requires zero coordinate conversions or 2D searches.

    s = struct('linIdx', uint32([]), 'isOn', false([]), 't', [], 'count', 0);
    duration = 0;

    if ~isfile(filePath)
        altPath = fullfile(fileparts(mfilename('fullpath')), filePath);
        if isfile(altPath)
            filePath = altPath;
        else
            warning('File not found for %s: %s', label, filePath);
            return;
        end
    end

    raw = readmatrix(filePath);
    assert(size(raw, 2) >= 4 && ~isempty(raw) && all(isfinite(raw(:, 1:4)), 'all'), ...
        'File "%s" must contain at least four finite numeric columns: x y polarity t_us.', filePath);

    x = round(raw(:, 1));
    y = round(raw(:, 2));
    if isZeroBased
        x = x + 1;
        y = y + 1;
    end

    pol = raw(:, 3);
    if all(pol == 0 | pol == 1)
        pol = 2 * pol - 1;
    elseif ~all(pol == -1 | pol == 1)
        error('Polarity column in "%s" must be either 0/1 or -1/+1.', filePath);
    end

    t_us = raw(:, 4);
    origin_us = min(t_us);
    duration = (max(t_us) - origin_us + 1) * 1e-6;
    t_s = (t_us - origin_us) * 1e-6;

    % Filter out-of-bounds pixels once at load time
    valid = (x >= 1 & x <= nx & y >= 1 & y <= ny);
    x = x(valid);
    y = y(valid);
    pol = pol(valid);
    t_s = t_s(valid);

    % Sort chronologically once
    [t_s, sortIdx] = sort(t_s);
    x = x(sortIdx);
    y = y(sortIdx);
    pol = pol(sortIdx);

    % Pre-compute 1D linear pixel index: y + (x - 1)*ny
    s.linIdx = uint32(y + (x - 1) * ny);
    s.isOn   = (pol > 0);
    s.t      = t_s;
    s.count  = numel(t_s);

    fprintf('Loaded & indexed %s (%s): %d valid events over %.6f s.\n', ...
        label, filePath, s.count, duration);
end

function play_realtime_comparison(streamA, durA, titleA, subA, ...
                                  streamB, durB, titleB, subB, ...
                                  nx, ny, window_s, playbackSpeed, targetFPS)

    totalDuration = max(durA, durB);
    window        = min(window_s, totalDuration);
    lastStart     = max(0, totalDuration - window);
    current       = 0;
    playing       = false;
    token         = 0;
    nPix          = nx * ny;

    % UI Colors (uint8 0-255 for fast direct RGB buffer writes)
    bgFig       = [0.96 0.97 0.96];
    bgCard      = [0.07 0.21 0.18];
    bgSensorU8  = uint8([10, 20, 18]);     % Near-black sensor background
    txtMain     = [1.00 1.00 1.00];
    txtSub      = [0.68 0.80 0.76];

    % Pre-build a 256x256 RGB Look-Up Table (LUT) for (onCount, offCount) -> uint8 RGB
    % Index 1 = 0 events, Index 256 = 255+ events
    countsVec = (0:255)';
    wVec      = min(1, log2(1 + countsVec) / 2.5);
    [wOffGrid, wOnGrid] = meshgrid(wVec, wVec); % rows: ON (1..256), cols: OFF (1..256)

    lutR = uint8(min(255, round(double(bgSensorU8(1)) + wOnGrid * 245 + wOffGrid * 71)));
    lutG = uint8(min(255, round(double(bgSensorU8(2)) + wOnGrid * 209 + wOffGrid * 204)));
    lutB = uint8(min(255, round(double(bgSensorU8(3)) + wOnGrid * 156 + wOffGrid * 166)));

    % Pre-allocate reusable frame buffers (avoids memory allocation during playback)
    rgbBufA = zeros(ny, nx, 3, 'uint8');
    rgbBufB = zeros(ny, nx, 3, 'uint8');
    rgbBufA(:,:,1) = bgSensorU8(1); rgbBufA(:,:,2) = bgSensorU8(2); rgbBufA(:,:,3) = bgSensorU8(3);
    rgbBufB(:,:,1) = bgSensorU8(1); rgbBufB(:,:,2) = bgSensorU8(2); rgbBufB(:,:,3) = bgSensorU8(3);

    fig = figure('Color', bgFig, 'Name', 'Real-Time Raw Event Stream Comparison', ...
        'NumberTitle', 'off', 'Position', [100 120 1380 680], ...
        'CloseRequestFcn', @close_view);

    % Top Banner Labels
    uicontrol(fig, 'Style', 'text', 'String', 'Event Stream Comparison', ...
        'Units', 'normalized', 'Position', [0.025 0.935 0.45 0.04], ...
        'BackgroundColor', bgFig, 'ForegroundColor', [0.15 0.20 0.18], ...
        'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'left');
    fpsLabel = uicontrol(fig, 'Style', 'text', ...
        'String', sprintf('Full sensor · %.1fx Real-Time', playbackSpeed), ...
        'Units', 'normalized', 'Position', [0.525 0.935 0.45 0.04], ...
        'BackgroundColor', bgFig, 'ForegroundColor', [0.40 0.48 0.45], ...
        'FontSize', 10, 'HorizontalAlignment', 'right');

    % Create Left and Right Display Cards
    [imA, footA] = create_stream_card([0.025 0.13 0.462 0.79], titleA, subA, rgbBufA);
    [imB, footB] = create_stream_card([0.513 0.13 0.462 0.79], titleB, subB, rgbBufB);

    % Bottom Playback Bar
    playBtn = uicontrol(fig, 'Style', 'pushbutton', 'String', 'Play', ...
        'Units', 'normalized', 'Position', [0.025 0.035 0.08 0.055], ...
        'FontSize', 10, 'FontWeight', 'bold', 'Callback', @play_pause);

    slider = uicontrol(fig, 'Style', 'slider', 'Units', 'normalized', ...
        'Position', [0.12 0.042 0.855 0.04], 'Min', 0, 'Max', max(lastStart, eps), ...
        'Value', 0, 'SliderStep', [min(1, 0.025/max(lastStart, eps)), ...
        min(1, 0.250/max(lastStart, eps))], 'Callback', @seek);
    if lastStart == 0, set(slider, 'Enable', 'off'); end

    draw_frame();

    % ---- Nested Helper Functions ----------------------------------------
    function [imHandle, footerHandle] = create_stream_card(pos, cardTitle, cardSub, initRGB)
        card = uipanel(fig, 'Units', 'normalized', 'Position', pos, ...
            'BackgroundColor', bgCard, 'BorderType', 'none');

        uicontrol(card, 'Style', 'text', 'String', cardTitle, ...
            'Units', 'normalized', 'Position', [0.035 0.91 0.93 0.055], ...
            'BackgroundColor', bgCard, 'ForegroundColor', txtMain, ...
            'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'left');
        uicontrol(card, 'Style', 'text', 'String', cardSub, ...
            'Units', 'normalized', 'Position', [0.035 0.865 0.93 0.045], ...
            'BackgroundColor', bgCard, 'ForegroundColor', txtSub, ...
            'FontSize', 9.5, 'HorizontalAlignment', 'left');

        ax = axes('Parent', card, 'Units', 'normalized', ...
            'Position', [0.005 0.085 0.99 0.765]);
        imHandle = image(ax, initRGB);
        axis(ax, 'image', 'off');
        set(ax, 'XLim', [0.5 nx+0.5], 'YLim', [0.5 ny+0.5], 'YDir', 'reverse');

        footerHandle = uicontrol(card, 'Style', 'text', 'String', '', ...
            'Units', 'normalized', 'Position', [0.035 0.018 0.93 0.045], ...
            'BackgroundColor', bgCard, 'ForegroundColor', txtSub, ...
            'FontSize', 9.5, 'HorizontalAlignment', 'left');
    end

    function [rgb, nRaw] = render_fast(s, coverage, t0, t1, rgb)
        % Reset background in-place
        rgb(:,:,1) = bgSensorU8(1);
        rgb(:,:,2) = bgSensorU8(2);
        rgb(:,:,3) = bgSensorU8(3);

        if s.count == 0 || t0 >= coverage
            nRaw = 0;
            return;
        end

        % O(log N) binary search for start and end event indices
        iStart = lower_bound(s.t, t0);
        iEnd   = lower_bound(s.t, min(t1, coverage)) - 1;

        if iEnd < iStart
            nRaw = 0;
            return;
        end

        nRaw = iEnd - iStart + 1;
        idxSlice = s.linIdx(iStart:iEnd);
        onSlice  = s.isOn(iStart:iEnd);

        % Fast 1D accumulation on only the active pixels
        onCounts  = min(255, accumarray(idxSlice(onSlice),  1, [nPix 1]));
        offCounts = min(255, accumarray(idxSlice(~onSlice), 1, [nPix 1]));

        active = find(onCounts > 0 | offCounts > 0);
        if isempty(active)
            return;
        end

        % Map active pixels through the precomputed 256x256 LUT
        lutIdx = onCounts(active) + 1 + offCounts(active) * 256;
        rgb(active)          = lutR(lutIdx);
        rgb(active + nPix)   = lutG(lutIdx);
        rgb(active + 2*nPix) = lutB(lutIdx);
    end

    function idx = lower_bound(vec, val)
        % Fast O(log N) binary search for first element >= val
        lo = 1;
        hi = numel(vec) + 1;
        while lo < hi
            mid = floor((lo + hi) * 0.5);
            if vec(mid) < val
                lo = mid + 1;
            else
                hi = mid;
            end
        end
        idx = lo;
    end

    function draw_frame()
        if ~isgraphics(fig), return; end
        finish = min(totalDuration, current + window);

        [rgbBufA, nA] = render_fast(streamA, durA, current, finish, rgbBufA);
        [rgbBufB, nB] = render_fast(streamB, durB, current, finish, rgbBufB);

        set(imA, 'CData', rgbBufA);
        set(imB, 'CData', rgbBufB);

        set(footA, 'String', sprintf('%.3f–%.3f s · %d raw events', current, finish, nA));
        set(footB, 'String', sprintf('%.3f–%.3f s · %d raw events', current, finish, nB));
        set(slider, 'Value', current);
    end

    function seek(source, ~)
        playing = false; token = token + 1; set(playBtn, 'String', 'Play');
        current = get(source, 'Value');
        draw_frame();
    end

    function play_pause(~, ~)
        if playing
            playing = false; token = token + 1; set(playBtn, 'String', 'Play');
            return;
        end
        if current >= lastStart, current = 0; end
        playing = true; token = token + 1; myToken = token; set(playBtn, 'String', 'Pause');

        % Anchor simulation time to true wall-clock time
        tStartSim  = current;
        tWallStart = tic;
        frameCount = 0;
        minFrameDt = 1 / targetFPS;

        while isgraphics(fig) && playing && myToken == token
            tFrameStart = tic;

            % Advance stream time strictly by elapsed real-world wall-clock time
            current = min(lastStart, tStartSim + toc(tWallStart) * playbackSpeed);

            draw_frame();
            drawnow limitrate;
            frameCount = frameCount + 1;

            if ~isgraphics(fig) || ~playing || myToken ~= token, return; end
            if current >= lastStart, break; end

            % Update live FPS display every 15 frames
            if mod(frameCount, 15) == 0
                actualFPS = frameCount / toc(tWallStart);
                set(fpsLabel, 'String', sprintf('Full sensor · %.1fx Real-Time (%.0f FPS)', ...
                    playbackSpeed, actualFPS));
            end

            % Yield only the remaining time needed to cap at targetFPS
            frameElapsed = toc(tFrameStart);
            if frameElapsed < minFrameDt
                pause(minFrameDt - frameElapsed);
            else
                pause(0.001); % Minimal yield so UI buttons/slider stay responsive
            end
        end

        if isgraphics(fig) && myToken == token
            playing = false;
            set(playBtn, 'String', 'Play');
        end
    end

    function close_view(~, ~)
        playing = false; token = token + 1;
        delete(fig);
    end
end