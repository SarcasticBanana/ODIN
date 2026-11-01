function cfg = init_matched_filter(nw, nx, ny, steadyMask, boxHalf, cenHalf, faBudget)
    % INIT_MATCHED_FILTER Generates the optical stencil and packages settings.
    
    % Call external function to generate the 2D Point Spread Function (PSF) 
    % matrix (e.g., 21x21 pixels). This models the physical shape of the target's light.
    g = matchedfilter_psf(nw);
    
    % Running a full 2D convolution on every tracker loop is computationally heavy.
    % Singular Value Decomposition (SVD) splits the 2D template (g) into two 
    % 1D vectors (gx and gy). Sliding these 1D vectors sequentially produces 
    % the exact same tracking result but executes exponentially faster.
    [U,Sv,V] = svd(g);
    gy = U(:,1)*sqrt(Sv(1,1));
    gx = V(:,1)*sqrt(Sv(1,1));
    
    % Ensure the vectors didn't flip signs during decomposition
    if sum(gy) < 0, gy = -gy; gx = -gx; end
    
    % Package all the tracking bounds, noise budgets, and optical stencils 
    % into a single configuration structure to easily pass into the tracker.
    cfg = struct('nx',nx,'ny',ny,'mask',steadyMask,'gy',gy,'gx',gx, ...
                 'boxHalf',boxHalf,'cenHalf',cenHalf,'fa',faBudget);
end