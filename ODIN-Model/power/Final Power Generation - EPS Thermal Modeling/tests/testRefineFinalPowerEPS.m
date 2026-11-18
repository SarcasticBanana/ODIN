function report=testRefineFinalPowerEPS(cfg,out)
% TESTREFINEFINALPOWEREPS Private-grid tests against native unchanged ODIN input.
% Run odinPath, add this package parent and tests folder, then:
%   cfg=config(); out=propagate(cfg,1); testRefineFinalPowerEPS(cfg,out)
% The supplied native ephemeris is not re-propagated or changed by this test.
% No arguments uses config()/propagate() already selected from the ODIN checkout.
if nargin==0, cfg=config(); end
if nargin<2, out=propagate(cfg,1); end
beforeCfg=cfg; beforeOut=out; beforePath=path; beforeFolder=pwd;
shared={'config','propagate','eom','shadow','eclipseTimes','sunVector'};
beforeResolution=cellfun(@which,shared,'UniformOutput',false);
count=0;
[refined,info]=finalPowerEPS.refineOrbit(out,cfg);
[present,index]=ismember(out.t,refined.t);
check(all(present),'Every original time must remain in the analysis grid.');
check(isequal(refined.r(:,index),out.r),'Original position knots changed.');
check(isequal(refined.v(:,index),out.v),'Original velocity knots changed.');
check(all(diff(refined.t)>0),'Refinement produced repeated or decreasing times.');
check(numel(refined.t)>=numel(out.t),'Refinement discarded source times.');
check(refined.t(1)==out.t(1) && refined.t(end)==out.t(end), ...
    'The helper extended or shortened the supplied ephemeris.');
check(all(refined.t>=out.t(1) & refined.t<=out.t(end)), ...
    'Refinement extrapolated outside the ephemeris.');
check(~info.endpointExtended,'Refinement claims it extended the endpoint.');
check(isequal(refined.hvec,cross(refined.r,refined.v,1)), ...
    'Angular momentum was not recomputed on the refined grid.');
check(all(isfinite(refined.illum)) && all(refined.illum>=0 & refined.illum<=1), ...
    'Refined illumination is invalid.');
contacts=[refined.eclipseTimes.umbraIn,refined.eclipseTimes.umbraOut, ...
    refined.eclipseTimes.penumbraIn,refined.eclipseTimes.penumbraOut];
check(all(ismember(contacts,refined.t)),'A detected contact is absent from the analysis grid.');

% A requested endpoint later than the data is reported, never invented.
short=out; short.requestedDuration_s=out.t(end)+123;
[shortRefined,shortInfo]=finalPowerEPS.refineOrbit(short,cfg);
check(shortInfo.requestedDurationKnown && ...
    abs(shortInfo.requestedDurationShortfall_s-123)<1e-9, ...
    'Missing endpoint duration was not reported correctly.');
check(shortRefined.t(end)==out.t(end) && ~shortInfo.endpointExtended, ...
    'Missing endpoint data was extrapolated.');

% Static metadata remains; unknown sampled diagnostics cannot retain stale lengths.
metadata=out; n=numel(out.t);
metadata.privateSampledDiagnostic=1:n;
metadata.privateSampledFlag=true(1,n);
metadata.orbitSolution=struct('placeholder','not used by this helper');
metadata.staticLabel='ODIN test metadata';
[clean,metadataInfo]=finalPowerEPS.refineOrbit(metadata,cfg);
omitted={'privateSampledDiagnostic','privateSampledFlag','orbitSolution'};
check(all(ismember(omitted,metadataInfo.omittedFields)), ...
    'Dropped sampled fields were not named in the provenance.');
check(~any(isfield(clean,omitted)),'Stale sampled metadata survived refinement.');
check(strcmp(clean.staticLabel,metadata.staticLabel),'Static text metadata changed.');

noShadowCfg=cfg; noShadowCfg.use.eclipse=false;
[noShadow,noShadowInfo]=finalPowerEPS.refineOrbit(out,noShadowCfg);
check(all(noShadow.illum==1) && ~noShadowInfo.shadowEnabled, ...
    'Eclipse-disabled configuration did not remain fully illuminated.');
check(isempty(noShadow.eclipseTimes.umbraIn) && isempty(noShadow.eclipseTimes.umbraOut) && ...
    isempty(noShadow.eclipseTimes.penumbraIn) && isempty(noShadow.eclipseTimes.penumbraOut), ...
    'Eclipse-disabled refinement generated contact events.');
check(isequal(noShadow.t,out.t),'Eclipse-disabled refinement changed the grid unnecessarily.');

reject(@()finalPowerEPS.refineOrbit(out,cfg,struct('unknownOption',1)), ...
    'finalPowerEPS:refineOrbit:option','Unknown helper options must be rejected.');
reject(@()finalPowerEPS.run(cfg,out,struct('unknownOption',1)), ...
    'finalPowerEPS:option','Unknown run options must be rejected.');

% Construct a cylindrical boundary crossing with a contact at the midpoint.
% This is a synthetic geometry check only, not a replacement orbit model.
% Centering it on the evaluated Sun direction makes the known contact exact.
contactCfg=cfg; contactCfg.use.eclipse=true; contactCfg.eclipse.model='cylindrical';
[~,sun]=finalPowerEPS.sunVector(1,contactCfg); sun=sun(:,1);
axis=cross(sun,[0;0;1]);
if norm(axis)<1e-8, axis=cross(sun,[0;1;0]); end
axis=axis/norm(axis);
center=-(contactCfg.Re+600)*sun+contactCfg.Re*axis;
contactOrbit=struct('t',[0 2],'r',[center-100*axis,center+100*axis], ...
    'v',[100*axis,100*axis],'T',2);
[contactRefined,~]=finalPowerEPS.refineOrbit(contactOrbit,contactCfg);
contactTimes=[contactRefined.eclipseTimes.umbraIn,contactRefined.eclipseTimes.umbraOut];
check(~isempty(contactTimes) && min(abs(contactTimes-1))<1e-8, ...
    'Synthetic cylindrical contact was not found at its known midpoint.');
check(numel(contactRefined.t)>numel(contactOrbit.t) && all(ismember(contactTimes,contactRefined.t)), ...
    'The known contact did not add samples.');
check(all(diff(contactRefined.t)>0),'One-sided cylindrical samples are not strictly increasing.');

% Opposite normals cannot be linearly normalized at the contact midpoint.
contactOrbit.panelNormal_ECI=[sun,-sun];
reject(@()finalPowerEPS.refineOrbit(contactOrbit,contactCfg), ...
    'finalPowerEPS:refineOrbit:normals','Antipodal attitude interpolation must be rejected.');

% Ordinary supplied normals retain their original exact knot values.
constantNormals=out; constantNormals.panelNormal_ECI=repmat([0;0;1],1,n);
[normalRefined,~]=finalPowerEPS.refineOrbit(constantNormals,cfg);
[~,normalIndex]=ismember(out.t,normalRefined.t);
check(isequal(normalRefined.panelNormal_ECI(:,normalIndex),constantNormals.panelNormal_ECI), ...
    'Supplied normal values changed at source times.');
check(max(abs(vecnorm(normalRefined.panelNormal_ECI)-1))<1e-12, ...
    'Refined panel normals lost unit magnitude.');

check(isequaln(cfg,beforeCfg) && isequaln(out,beforeOut),'Refinement changed native inputs.');
check(strcmp(path,beforePath) && strcmp(pwd,beforeFolder),'Refinement changed path or current folder.');
check(isequal(beforeResolution,cellfun(@which,shared,'UniformOutput',false)), ...
    'Refinement shadowed a shared team function.');
report=struct('checksPassed',count,'sourceSamples',numel(out.t), ...
    'analysisSamples',numel(refined.t),'sourceEndpoint_s',out.t(end), ...
    'analysisEndpoint_s',refined.t(end),'detectedContacts',numel(contacts));
fprintf('PASS: %d private-grid, ephemeris-preservation and option/attitude checks.\n',count);
    function check(ok,message)
        assert(all(ok(:)),message); count=count+1;
    end
    function reject(action,identifier,message)
        rejected=false;
        try
            action();
        catch err
            rejected=strcmp(err.identifier,identifier);
        end
        check(rejected,message);
    end
end
