function steadyMask = build_steady_source_mask(t, x, y, T, nx, ny, nSeg, segFrac, maskPad)
    % BUILD_STEADY_SOURCE_MASK Identifies and removes stationary noise.
    % Defective "hot pixels" fire continuously (100+ events/sec), and stationary 
    % stars create constant streaks. Without removal, the matched filter will 
    % lock onto these instead of faint orbital debris.
    
    % Chop the total recording duration (T) into a set number of equal time 
    % chunks (nSeg), assigning every event to a specific chunk.
    seg  = min(floor(t/T*nSeg)+1, nSeg);
    
    % Create a blank 2D grid matching the sensor resolution to tally activity
    nOcc = zeros(ny,nx);               
    
    % Loop through each time chunk. If a pixel fired at all during that chunk, 
    % add 1 to its lifetime tally.
    for s = 1:nSeg
        ii = (seg == s);
        if ~any(ii), continue; end
        nOcc = nOcc + (accumarray([y(ii) x(ii)],1,[ny nx]) > 0);
    end
    
    % If a pixel fired in a majority of the chunks (e.g., 6 out of 10), it 
    % is exhibiting steady-state noise behavior. Flag it.
    steady = nOcc >= ceil(segFrac*nSeg);
    
    % Dilate the flagged pixels outward by maskPad (e.g., a 2-pixel border) 
    % to ensure the edges of a bright star or hot pixel don't leak through.
    steadyMask = conv2(double(steady), ones(2*maskPad+1), 'same') > 0;
    fprintf('Steady-source pixels masked: %d\n', nnz(steadyMask));
end