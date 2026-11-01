function [px,py,keep,res] = robust_line_fit(tt,xx,yy,clipSigma,minResid)
% This fits a line to the extracted X,Y centroid coordinates over time to find velocity.
% To prevent a single bad data point from ruining the line, it uses a "Theil-Sen estimator",
% which calculates the slope between every single pair of points and takes the median slope.

    n   = numel(tt);
    sub = round(linspace(1,n,min(n,600)));      % Cap the max number of points to keep it fast
    ts = tt(sub);  xsub = xx(sub);  ysub = yy(sub);
    
    % Find the slope between all possible pairs of points
    [I,J] = find(triu(true(numel(sub)),1));
    dtp  = ts(J) - ts(I);
    good = dtp > 0;
    
    % The median slope is highly resistant to massive outliers
    vx0 = median((xsub(J(good)) - xsub(I(good)))./dtp(good));
    vy0 = median((ysub(J(good)) - ysub(I(good)))./dtp(good));
    px = [vx0, median(xx - vx0*tt)];
    py = [vy0, median(yy - vy0*tt)];

    % Iterate a standard Least Squares fit, but aggressively delete points (Sigma Clipping) 
    % that fall too far away from the median line we just calculated.
    keep = true(size(tt));
    for it = 1:20
        res = hypot(xx - polyval(px,tt), yy - polyval(py,tt)); % Distance from point to line
        s = max(1.4826*median(res(keep)), minResid);
        newKeep = res < clipSigma*s; % Flag points outside the 3-sigma bound
        if it > 1 && isequal(newKeep,keep), break; end
        keep = newKeep;
        px = polyfit(tt(keep),xx(keep),1); % Re-fit line using only the good points
        py = polyfit(tt(keep),yy(keep),1);
    end
    res = hypot(xx - polyval(px,tt), yy - polyval(py,tt));
end