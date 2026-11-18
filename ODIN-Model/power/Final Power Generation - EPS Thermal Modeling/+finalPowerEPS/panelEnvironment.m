function e=panelEnvironment(out,cfg,p,angles_deg)
% PANELENVIRONMENT Face-resolved Earth IR/albedo from diffuse-sphere geometry.
% Fluxes are absorbed W/m2 of ONE physical panel face. IR is not electrical input.
% Quadrature integrates radiance over Earth's apparent disk at the spacecraft:
% F=(1/pi)*integral(max(n dot ray,0) dOmega), L_IR=earthInfrared/pi,
% L_albedo=albedo*S_unattenuated*max(nEarth dot Sun,0)/pi.
% This includes limb projection, visible illuminated terrain, and inverse-square
% geometry; it does not multiply albedo by the satellite's eclipse fraction.
% out.panelNormal_ECI may supply measured/attitude-derived normals (3 x N).
% Otherwise use Sun tracking, with fixed-angle tilt toward the orbit normal.
angles_deg=angles_deg(:); n=numel(out.t); m=numel(angles_deg);
names={'frontEarthViewFactor','backEarthViewFactor', ...
    'frontAbsorbedEarthFlux_W_m2','backAbsorbedEarthFlux_W_m2', ...
    'frontEarthIR_W_m2','backEarthIR_W_m2','frontAlbedo_W_m2','backAlbedo_W_m2'};
for k=1:numel(names), e.(names{k})=zeros(m,n); end
e.enabled=p.environment.enabled;
e.settings=p.environment;
e.normalSource='Sun tracking; fixed-angle tilt toward projected orbit normal';
e.structurePathsEnabled=p.mountConductancePerArea_W_m2K>0 || ...
    p.frontStructureViewFactor>0 || p.backStructureViewFactor>0;
if ~e.enabled, e.normalSource='Environment disabled explicitly'; return; end
if ~isfield(out,'r') || ~isfield(out,'rsun')
    error('panelEnvironment:geometry','Earth environment requires out.r and out.rsun; explicitly disable environment for geometry-free tests.');
end
validateattributes(out.r,{'numeric'},{'real','finite','size',[3 n]});
validateattributes(out.rsun,{'numeric'},{'real','finite','size',[3 n]});
v=p.environment;
validateattributes(v.earthInfrared_W_m2,{'numeric'},{'real','finite','scalar','nonnegative'});
validateattributes(v.earthAlbedo,{'numeric'},{'real','finite','scalar','>=',0,'<=',1});
validateattributes(v.radialOrder,{'numeric'},{'scalar','integer','>=',2});
validateattributes(v.azimuthCount,{'numeric'},{'scalar','integer','>=',8});
if ~strcmp(v.normalMode,'sun_tracking'), error('panelEnvironment:normal','Use sun_tracking or provide out.panelNormal_ECI.'); end
supplied=isfield(out,'panelNormal_ECI');
if supplied
    validateattributes(out.panelNormal_ECI,{'numeric'},{'real','finite','size',[3 n]});
    if any(abs(vecnorm(out.panelNormal_ECI)-1)>1e-6)
        error('panelEnvironment:normal','Panel normals must be unit vectors.');
    end
    if m~=1, error('panelEnvironment:normal','Supplied attitude normals describe a single case, not an independent fixed-angle sweep.'); end
    e.normalSource='Supplied out.panelNormal_ECI';
end
[nodes,weights]=gaussRule(v.radialOrder);
phi=(0:v.azimuthCount-1)*2*pi/v.azimuthCount;
for k=1:n
    r=out.r(:,k); radius=norm(r);
    if radius<=cfg.Re, error('panelEnvironment:altitude','Spacecraft must be outside Earth.'); end
    z=-r/radius; shat=out.rsun(:,k)/norm(out.rsun(:,k));
    axis=cross(z,shat);
    if norm(axis)<1e-10, axis=cross(z,[0;0;1]); end
    if norm(axis)<1e-10, axis=cross(z,[0;1;0]); end
    x=axis/norm(axis); y=cross(z,x);
    muMin=sqrt(1-(cfg.Re/radius)^2);
    mu=muMin+(nodes+1)*(1-muMin)/2;
    w=weights*(1-muMin)/v.azimuthCount; % dOmega/pi = dmu*dphi/pi
    [MU,PHI]=ndgrid(mu,phi);
    rays=z*MU(:).'+x*(sqrt(1-MU(:).^2).*cos(PHI(:))).'+ ...
        y*(sqrt(1-MU(:).^2).*sin(PHI(:))).';
    omega=repmat(w,1,v.azimuthCount); omega=omega(:).';
    distance=radius*MU(:).'-sqrt(max(0,cfg.Re^2-radius^2*(1-MU(:).'.^2)));
    terrain=(r+rays.*distance)/cfg.Re;
    sunCos=max(0,shat.'*terrain);
    if supplied
        normals=out.panelNormal_ECI(:,k).';
    else
        if isfield(out,'v'), tilt=cross(r,out.v(:,k)); else, tilt=[0;0;1]; end
        tilt=tilt-shat*dot(tilt,shat);
        if norm(tilt)<1e-10, tilt=cross(shat,[0;1;0]); end
        if norm(tilt)<1e-10, tilt=cross(shat,[1;0;0]); end
        tilt=tilt/norm(tilt);
        normals=cosd(angles_deg)*shat.'+sind(angles_deg)*tilt.';
    end
    front=max(0,normals*rays); back=max(0,-normals*rays);
    fFront=front*omega.'; fBack=back*omega.';
    % Existing structure view factors mask that fraction of the visible sky.
    fFront=fFront*(1-p.frontStructureViewFactor);
    fBack=fBack*(1-p.backStructureViewFactor);
    e.frontEarthViewFactor(:,k)=fFront;
    e.backEarthViewFactor(:,k)=fBack;
    e.frontEarthIR_W_m2(:,k)=p.frontEmissivity*v.earthInfrared_W_m2*fFront;
    e.backEarthIR_W_m2(:,k)=p.backEmissivity*v.earthInfrared_W_m2*fBack;
    solar=cfg.S0*(cfg.AU/norm(out.rsun(:,k)))^2;
    e.frontAlbedo_W_m2(:,k)=p.solarAbsorptivity*v.earthAlbedo*solar* ...
        (front*(omega.*sunCos).')*(1-p.frontStructureViewFactor);
    e.backAlbedo_W_m2(:,k)=p.solarAbsorptivity*v.earthAlbedo*solar* ...
        (back*(omega.*sunCos).')*(1-p.backStructureViewFactor);
end
e.frontAbsorbedEarthFlux_W_m2=e.frontEarthIR_W_m2+e.frontAlbedo_W_m2;
e.backAbsorbedEarthFlux_W_m2=e.backEarthIR_W_m2+e.backAlbedo_W_m2;
e.assumptions=['Spherical Earth, uniform 237 W/m2 IR and 0.30 Lambertian albedo by default; ' ...
    'no weather/map model. Structure view fractions uniformly mask the sky, not resolved occlusion. ' ...
    'Structure temperatures are prescribed; unknown mounting conductance remains zero.'];
end
function [x,w]=gaussRule(n)
% Gauss-Legendre weights on [-1,1], cached by quadrature order.
persistent lastN lastX lastW
if isequal(lastN,n), x=lastX; w=lastW; return; end
b=(1:n-1)./sqrt(4*(1:n-1).^2-1);
[vectors,values]=eig(diag(b,1)+diag(b,-1));
[x,index]=sort(diag(values)); w=2*vectors(1,index).'.^2;
lastN=n; lastX=x; lastW=w;
end
