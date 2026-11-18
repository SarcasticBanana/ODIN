function q=panelEnvironmentAt(p,e,t,queryTime,method)
% PANELENVIRONMENTAT Apply precomputed face boundary histories to kernel settings.
q=p;
if ~e.enabled, return; end
fields={'frontEarthViewFactor','backEarthViewFactor', ...
    'frontAbsorbedEarthFlux_W_m2','backAbsorbedEarthFlux_W_m2'};
% ODE calls request one instant repeatedly. Rebuilding/transposing the entire
% angle-by-time matrix in interp1 at every stage is unnecessary: locate the
% interval once and combine only its two columns for all faces and angles.
fastScalar=isscalar(queryTime) && numel(t)>1 && ...
    any(strcmp(method,{'linear','previous','next','nearest'}));
if fastScalar
    lo=1; hi=numel(t)-1;
    while lo<hi
        middle=floor((lo+hi+1)/2);
        if queryTime>=t(middle), lo=middle; else, hi=middle-1; end
    end
    weight=(queryTime-t(lo))/(t(lo+1)-t(lo));
    if strcmp(method,'previous'), weight=double(queryTime>=t(lo+1)); end
    if strcmp(method,'next'), weight=double(queryTime>t(lo)); end
    if strcmp(method,'nearest'), weight=double(weight>=0.5); end
end
for k=1:numel(fields)
    name=fields{k};
    if isequal(t,queryTime)
        q.(name)=e.(name);
    elseif fastScalar
        q.(name)=(1-weight)*e.(name)(:,lo)+weight*e.(name)(:,lo+1);
    else
        q.(name)=interp1(t,e.(name).',queryTime,method).';
    end
end
end
