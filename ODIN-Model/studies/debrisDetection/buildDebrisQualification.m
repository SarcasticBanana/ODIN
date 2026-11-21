function qualify = buildDebrisQualification(masterSizeFile, sensorTable)
% Build once before the orbit loop. sensorTable is the corrected readtable output.
% P = qualify(altitudeKm,rangeIndex,speedKmS) returns the qualifying fraction
% of the local 1-5 cm population. Inputs are equal-sized or scalar.
% rangeIndex indexes unique(sensorTable.SlantRange_km) in ascending order.
% speedKmS is apparent transverse speed in the camera frame.
% Diameter fractions are assumed independent of inclination and declination
% within each 10 km altitude bin. Outside the modeled shell P is zero.
    data = readmatrix(masterSizeFile,'FileType','text','CommentStyle','#');
    data = data(isfinite(data(:,1)) & data(:,2)>=0.01 & data(:,2)<0.05,:);
    [altitude,~,ih] = unique(data(:,1));
    [diameter,~,id] = unique(100*data(:,2)); % meters to centimeters
    nH = numel(altitude);
    nD = numel(diameter);
    assert(all(accumarray([ih id],1,[nH nD])==1,'all'), ...
        'MASTER must contain one row per altitude and diameter.');
    density = accumarray([ih id],data(:,end),[nH nD]);
    assert(all(isfinite(density) & density>=0,'all') && all(sum(density,2)>0), ...
        'Invalid MASTER density or empty altitude bin.');
    % Column j gives the fraction at or above the j-th bin center.
    tail = [fliplr(cumsum(fliplr(density),2))./sum(density,2), zeros(nH,1)];
    assert(all(diff(altitude)==10),'Expected 10 km MASTER altitude bins.');
    edges = [altitude-5; altitude(end)+5];

    [ranges,~,ir] = unique(sensorTable.SlantRange_km);
    [sensorD,~,is] = unique(sensorTable.DebrisSize_cm);
    assert(all(accumarray([ir is],1)==1,'all'),'Expected one sensor configuration.');
    limit = accumarray([ir is],sensorTable.MaxRelativeVel_kms);
    % Interpolate velocity limits across diameter at each supplied range.
    limit = interp1(sensorD,limit.',diameter,'linear').';
    assert(all(isfinite(limit),'all') && all(diff(limit,1,2)>=0,'all'), ...
        'Sensor thresholds must span MASTER diameters and increase with size.');
    qualify = @evaluate;

    function p = evaluate(h,k,v)
        assert(isequal(size(h),size(k)) || isscalar(h) || isscalar(k));
        shape = size(h + zeros(size(k)));
        assert(isequal(shape,size(v)) || prod(shape)==1 || isscalar(v));
        template = zeros(size(zeros(shape)+v));
        h = h + template;
        k = k + template;
        v = v + template;
        assert(all(k>=1 & k<=numel(ranges) & k==fix(k),'all'));
        assert(all(isfinite(v) & v>=0,'all'),'Expected finite nonnegative speed.');
        bin = discretize(h,edges);
        active = isfinite(bin);
        bin = bin(active);
        r = k(active);
        speed = v(active);
        % Binary search for the first qualifying diameter, inclusive at equality.
        lo = ones(size(speed));
        hi = (nD+1)*lo;
        while any(lo<hi,'all')
            moving = lo<hi;
            mid = floor((lo+hi)/2);
            fast = limit(sub2ind(size(limit),r,min(mid,nD))) < speed;
            lo(moving & fast) = mid(moving & fast)+1;
            hi(moving & ~fast) = mid(moving & ~fast);
        end
        p = template;
        p(active) = tail(sub2ind(size(tail),bin,lo));
    end
end
