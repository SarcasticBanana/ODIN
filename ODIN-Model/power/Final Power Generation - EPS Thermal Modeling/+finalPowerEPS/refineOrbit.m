function [refined,info]=refineOrbit(out,cfg,options)
% REFINEORBIT Private power-analysis grid from an unchanged team orbit.
% Positions/velocities use cubic Hermite interpolation of supplied km/km/s
% states. No equations of motion are substituted and no endpoint is extended.
% Original r/v samples and supplied panel normals are retained exactly.
% Sun and shadow are evaluated independently on the analysis grid. Unknown
% numeric time-series fields are omitted and named in info.omittedFields;
% static metadata is retained. No path, base workspace, or input is modified.
if nargin<3, options=struct(); end
validateattributes(out,{'struct'},{'scalar'});
validateattributes(cfg,{'struct'},{'scalar'});
validateattributes(options,{'struct'},{'scalar'});
allowedOptions={'searchStep_s','transitionStep_s'};
for optionName=fieldnames(options).'
    if ~ismember(optionName{1},allowedOptions)
        error('finalPowerEPS:refineOrbit:option','Unknown refinement option: %s.',optionName{1});
    end
end
for field={'Re','Rsun','AU','epoch'}
    if ~isfield(cfg,field{1})
        error('finalPowerEPS:refineOrbit:config','Configuration requires %s.',field{1});
    end
end
searchStep=10; transitionStep=1;
if isfield(options,'searchStep_s'), searchStep=options.searchStep_s; end
if isfield(options,'transitionStep_s'), transitionStep=options.transitionStep_s; end
validateattributes(searchStep,{'numeric'},{'real','finite','scalar','positive'});
validateattributes(transitionStep,{'numeric'},{'real','finite','scalar','positive'});
for field={'t','r','v'}
    if ~isfield(out,field{1}), error('finalPowerEPS:refineOrbit:input','Orbit output requires %s.',field{1}); end
end
t=double(out.t(:).');
validateattributes(t,{'numeric'},{'real','finite','vector','nonempty'});
if numel(t)<2 || any(diff(t)<=0)
    error('finalPowerEPS:refineOrbit:time','Supply at least two strictly increasing orbit times.');
end
n=numel(t);
validateattributes(out.r,{'numeric'},{'real','finite','size',[3 n]});
validateattributes(out.v,{'numeric'},{'real','finite','size',[3 n]});
validateattributes(cfg.Re,{'numeric'},{'real','finite','scalar','positive'});
if any(vecnorm(out.r)<=cfg.Re)
    error('finalPowerEPS:refineOrbit:altitude','All source positions must lie outside Earth.');
end
if isfield(out,'cfg') && isstruct(out.cfg) && isfield(out.cfg,'epoch') && ...
        ~isequal(out.cfg.epoch,cfg.epoch)
    error('finalPowerEPS:refineOrbit:epoch','cfg.epoch must match the supplied orbit epoch.');
end
enabled=true;
if isfield(cfg,'use') && isfield(cfg.use,'eclipse'), enabled=logical(cfg.use.eclipse); end
model='conical';
if isfield(cfg,'eclipse') && isfield(cfg.eclipse,'model'), model=char(cfg.eclipse.model); end
if ~any(strcmpi(model,{'conical','cylindrical'}))
    error('finalPowerEPS:refineOrbit:model','Unsupported shadow model.');
end
shadowCfg=cfg; shadowCfg.eclipse.model=model;

survey=unique([t,t(1):searchStep:t(end),t(end)]);
[rSurvey,~]=stateAt(survey);
sunSurvey=finalPowerEPS.sunVector(survey,cfg);
events=emptyEvents();
if enabled
    [gu,gp]=margins(rSurvey,sunSurvey);
    events.umbraIn=contacts(gu,1,-1);
    events.umbraOut=contacts(gu,1,1);
    events.penumbraIn=contacts(gp,2,-1);
    events.penumbraOut=contacts(gp,2,1);
    events.umbraDur=durations(events.umbraIn,events.umbraOut);
    events.penumbraDur=durations(events.penumbraIn,events.penumbraOut);
end
allContacts=unique([events.umbraIn,events.umbraOut,events.penumbraIn,events.penumbraOut]);
pieces=unique([t(1),allContacts,t(end)]);
transition=[];
if enabled && strcmpi(model,'conical')
    for k=1:numel(pieces)-1
        a=pieces(k); b=pieces(k+1); mid=(a+b)/2;
        [rm,~]=stateAt(mid);
        fraction=finalPowerEPS.shadow(rm,finalPowerEPS.sunVector(mid,cfg),shadowCfg);
        if fraction>0 && fraction<1
            transition=[transition,a:transitionStep:b,b]; %#ok<AGROW>
        end
    end
elseif enabled
    % A binary eclipse has a jump: retain one-sided samples on both sides.
    % This does not replace the shared shadow implementation or source output.
    epsilon=min(1e-3,transitionStep/100);
    transition=[allContacts-epsilon,allContacts+epsilon];
    transition=transition(transition>=t(1) & transition<=t(end));
end
analysisTime=unique([t,allContacts,transition]);
[r,v]=stateAt(analysisTime);
if any(vecnorm(r)<=cfg.Re)
    error('finalPowerEPS:refineOrbit:interpolation','Interpolated states enter Earth; provide a finer source ephemeris.');
end

refined=out;
known={'t','r','v','hvec','rsun','sunhat','illum','panelNormal_ECI','eclipseTimes'};
omitted={};
fields=fieldnames(out);
for k=1:numel(fields)
    name=fields{k}; value=out.(name);
    sampled=(isnumeric(value)||islogical(value)) && ~isscalar(value) && ...
        (size(value,2)==n || (isvector(value) && numel(value)==n));
    if (~ismember(name,known) && sampled) || strcmp(name,'orbitSolution')
        refined=rmfield(refined,name); omitted{end+1}=name; %#ok<AGROW>
    end
end
refined.t=analysisTime; refined.r=r; refined.v=v;
refined.hvec=cross(r,v,1);
[refined.rsun,refined.sunhat]=finalPowerEPS.sunVector(analysisTime,cfg);
if enabled
    refined.illum=finalPowerEPS.shadow(r,refined.rsun,shadowCfg);
else
    refined.illum=ones(size(analysisTime));
end
refined.eclipseTimes=events;
if isfield(out,'panelNormal_ECI')
    validateattributes(out.panelNormal_ECI,{'numeric'},{'real','finite','size',[3 n]});
    if any(abs(vecnorm(out.panelNormal_ECI)-1)>1e-6)
        error('finalPowerEPS:refineOrbit:normals','Supplied panel normals must be unit vectors.');
    end
    normals=interp1(t,out.panelNormal_ECI.',analysisTime,'linear').';
    magnitude=vecnorm(normals);
    if any(magnitude<1e-8)
        error('finalPowerEPS:refineOrbit:normals','Interpolated attitude is ambiguous; provide finer panel-normal samples.');
    end
    normals=normals./magnitude;
    [~,index]=ismember(t,analysisTime); normals(:,index)=out.panelNormal_ECI;
    refined.panelNormal_ECI=normals;
end
info=struct('enabled',true,'shadowEnabled',enabled,'sourceSamples',n,'analysisSamples',numel(analysisTime), ...
    'sourceStart_s',t(1),'sourceEnd_s',t(end),'analysisEnd_s',analysisTime(end), ...
    'endpointExtended',false,'searchStep_s',searchStep,'transitionStep_s',transitionStep, ...
    'interpolation','Cubic Hermite using supplied positions and velocities', ...
    'omittedFields',{omitted},'requestedDurationKnown',false,'requestedDurationShortfall_s',NaN);
if isfield(out,'requestedDuration_s')
    validateattributes(out.requestedDuration_s,{'numeric'},{'real','finite','scalar','nonnegative'});
    % Source times are elapsed seconds from the original simulation epoch.
    % requestedDuration_s is its requested elapsed endpoint, as in propagate.
    info.requestedDurationKnown=true;
    info.requestedDurationShortfall_s=max(0,out.requestedDuration_s-t(end));
end
info.limitations=['Refinement resolves the supplied ephemeris; it does not improve its dynamical accuracy. ' ...
    'Grazing events wholly missed by the search grid need finer source/search samples. ' ...
    'The supplied final time is retained; no missing orbital endpoint is extrapolated.'];

    function [position,velocity]=stateAt(query)
        query=query(:).';
        if any(query<t(1) | query>t(end))
            error('finalPowerEPS:refineOrbit:extrapolation','Query outside the supplied ephemeris.');
        end
        index=discretize(query,t); index(query==t(end))=n-1;
        h=t(index+1)-t(index); u=(query-t(index))./h;
        p0=out.r(:,index); p1=out.r(:,index+1);
        v0=out.v(:,index); v1=out.v(:,index+1);
        position=p0.*(2*u.^3-3*u.^2+1)+v0.*(h.*(u.^3-2*u.^2+u))+ ...
            p1.*(-2*u.^3+3*u.^2)+v1.*(h.*(u.^3-u.^2));
        velocity=p0.*((6*u.^2-6*u)./h)+v0.*(3*u.^2-4*u+1)+ ...
            p1.*((-6*u.^2+6*u)./h)+v1.*(3*u.^2-2*u);
        % Avoid roundoff changes at supplied knots.
        [exact,sourceIndex]=ismember(query,t);
        position(:,exact)=out.r(:,sourceIndex(exact));
        velocity(:,exact)=out.v(:,sourceIndex(exact));
    end
    function roots=contacts(g,kind,direction)
        left=g(1:end-1); right=g(2:end);
        if direction<0, crossing=left>=0 & right<0; else, crossing=left<0 & right>=0; end
        brackets=find(crossing); roots=zeros(1,numel(brackets));
        for j=1:numel(brackets)
            bracketIndex=brackets(j);
            roots(j)=fzero(@contactMargin,[survey(bracketIndex),survey(bracketIndex+1)]);
        end
        roots=unique(roots);
        function value=contactMargin(time)
            [position,~]=stateAt(time);
            [uMargin,pMargin]=margins(position,finalPowerEPS.sunVector(time,cfg));
            if kind==1, value=uMargin; else, value=pMargin; end
        end
    end
    function [uMargin,pMargin]=margins(position,sun)
        if strcmpi(model,'cylindrical')
            direction=sun./vecnorm(sun); along=sum(position.*direction,1);
            uMargin=max(along,vecnorm(position-along.*direction)-cfg.Re); pMargin=uMargin;
        else
            [separation,aSun,aEarth]=finalPowerEPS.shadowGeometry(position,sun,cfg);
            uMargin=separation-(aEarth-aSun); pMargin=separation-(aEarth+aSun);
        end
    end
end
function events=emptyEvents()
events=struct('umbraIn',[],'umbraOut',[],'penumbraIn',[],'penumbraOut',[], ...
    'umbraDur',[],'penumbraDur',[]);
end
function duration=durations(entry,exit)
duration=nan(size(entry));
for k=1:numel(entry)
    next=find(exit>entry(k),1);
    if ~isempty(next), duration(k)=exit(next)-entry(k); end
end
end
