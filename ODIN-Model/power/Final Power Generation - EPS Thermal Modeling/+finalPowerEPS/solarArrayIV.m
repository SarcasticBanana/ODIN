function result=solarArrayIV(inputs,overrides)
% SOLARARRAYIV 64-cell NanoAvionics DESIGN proxy using Spectrolab XTJ Prime cell behavior.
% This is NOT a vendor-verified NanoAvionics I-V curve. Ns/Np are assumed wiring.
% Default Ns=16, Np=4 uses all 64 cells; no 28 V voltage scale is imposed.
% Default mode enforces the existing panel available-power budget by current scaling.
% xtj_reference mode instead predicts the unmodified Spectrolab reference assembly.
% Geometry and coverage share arrayGeometry.m with the power/thermal models.

% Shared defaults and assumptions: edit finalPowerEPS.modelSettings('iv').
p=finalPowerEPS.modelSettings('iv');

if nargin==0, result.params=p; return; end
if nargin<2, overrides=struct(); end
if ~isstruct(inputs)||~isscalar(inputs)||~isstruct(overrides)||~isscalar(overrides)
    error('solarArrayIV:input','Inputs and overrides must be scalar structs.');
end
for name=fieldnames(overrides).'
    if ~isfield(p,name{1}), error('solarArrayIV:parameter','Unknown setting: %s',name{1}); end
    p.(name{1})=overrides.(name{1});
end
p.mode=validatestring(p.mode,{'power_consistent_proxy','xtj_reference'});
if isempty(p.parallelStrings), p.parallelStrings=p.totalCells/p.seriesCells; end
for name={'totalCells','seriesCells','parallelStrings', ...
        'referenceIrradiance_W_m2','cellVoc_V','cellVmp_V','cellJsc_A_cm2','cellJmp_A_cm2'}
    validateattributes(p.(name{1}),{'numeric'},{'real','finite','scalar','positive'});
end
if ~isempty(p.testVoltage_V), validateattributes(p.testVoltage_V,{'numeric'},{'real','finite','scalar','nonnegative'}); end
for name={'totalCells','seriesCells','parallelStrings','curvePoints','snapshotCount'}
    validateattributes(p.(name{1}),{'numeric'},{'real','finite','scalar','integer','positive'});
end
if p.seriesCells*p.parallelStrings~=p.totalCells
    error('solarArrayIV:topology','seriesCells*parallelStrings must equal totalCells.');
end
geometry=finalPowerEPS.arrayGeometry(p.arrayArea_m2,p.cellCoverage);
p.cellArea_cm2=geometry.totalCellArea_m2*1e4/p.totalCells;
validateattributes(p.coefficientRange_C,{'numeric'},{'real','finite','vector','numel',2});
if p.coefficientRange_C(2)<=p.coefficientRange_C(1), error('solarArrayIV:parameter','Invalid coefficient range.'); end
if p.curvePoints<3 || p.cellVmp_V>=p.cellVoc_V || p.cellJmp_A_cm2>=p.cellJsc_A_cm2
    error('solarArrayIV:parameter','Require at least 3 curve points, Vmp<Voc and Jmp<Jsc.');
end
for name={'referenceTemperature_C','cellVocTemperature_V_K','cellJscTemperature_A_cm2K', ...
        'cellVmpTemperature_V_K','cellJmpTemperature_A_cm2K'}
    validateattributes(p.(name{1}),{'numeric'},{'real','finite','scalar'});
end
validateattributes(p.referenceTemperature_C,{'numeric'},{'>',-273.15});
for name={'makePlots','showSummary'}
    validateattributes(double(p.(name{1})),{'numeric'},{'scalar','integer','>=',0,'<=',1});
end
validateattributes(p.figureSize_px,{'numeric'},{'real','finite','positive','vector','numel',2});
for name={'plotLineWidth','plotFontSize','plotPaddingFraction'}
    validateattributes(p.(name{1}),{'numeric'},{'real','finite','positive','scalar'});
end
allowed={'time_s','irradiance_W_m2','cellTemperature_C','availablePower_W'};
if any(~ismember(fieldnames(inputs),allowed)), error('solarArrayIV:input','Unknown input field.'); end
for name={'time_s','irradiance_W_m2','cellTemperature_C'}
    if ~isfield(inputs,name{1}), error('solarArrayIV:input','Missing %s.',name{1}); end
end
t=inputs.time_s(:).'; validateattributes(t,{'numeric'},{'real','finite','vector','nonempty'});
if any(diff(t)<=0), error('solarArrayIV:time','Times must increase.'); end
G=signal(inputs.irradiance_W_m2,numel(t)); T=signal(inputs.cellTemperature_C,numel(t));
if any(G<0|T<=-273.15), error('solarArrayIV:physics','Require nonnegative irradiance and temperatures above absolute zero.'); end
[aRef,rsRef,ilRef,i0Ref,fitResidual]=fitReference(p);
a=aRef*(T+273.15)/(p.referenceTemperature_C+273.15);
ilTemperature=ilRef+p.cellJscTemperature_A_cm2K*(T-p.referenceTemperature_C);
vocTemperature=p.cellVoc_V+p.cellVocTemperature_V_K*(T-p.referenceTemperature_C);
if any(ilTemperature<=0|vocTemperature<=0), error('solarArrayIV:temperature','Temperature extrapolation produces invalid cell parameters.'); end
% Rsh=infinity is the explicit reduced-order assumption. Effective diode thermal
% voltage a scales with absolute temperature; I0(T) reproduces the specified Voc(T).
x=vocTemperature./a;
logI0=log(ilTemperature)-x-log1p(-exp(-x));
i0=exp(logI0); il=ilTemperature.*G/p.referenceIrradiance_W_m2;
active=G>0;
voc=zeros(size(G)); voc(active)=a(active).*log1p(il(active)./i0(active));
if any(~isfinite(voc)), error('solarArrayIV:numerics','Unsupported extreme temperature/irradiance.'); end
lo=zeros(size(G)); hi=voc;
for k=1:60
    u=(lo+hi)/2; j=il-i0.*expm1(u./a); d=i0.*exp(u./a)./a;
    slope=j-d.*(u-2*j*rsRef);
    up=slope>0; lo(up)=u(up); hi(~up)=u(~up);
end
u=(lo+hi)/2; jm=max(0,il-i0.*expm1(u./a)); vm=max(0,u-jm*rsRef);
jsc=currentAt(zeros(size(G)),voc,il,i0,a,rsRef);
voltageScale=p.seriesCells;
currentScale=p.parallelStrings*p.cellArea_cm2*ones(size(G));
referencePmp=voltageScale*vm.*currentScale.*jm;
normalization=ones(size(G));
if strcmp(p.mode,'power_consistent_proxy')
    if ~isfield(inputs,'availablePower_W'), error('solarArrayIV:power','Proxy mode needs the shared panel availablePower_W history.'); end
    power=signal(inputs.availablePower_W,numel(t));
    if any(power<0 | (~active & power>1e-9)), error('solarArrayIV:power','Power must be nonnegative and zero in eclipse.'); end
    normalization=zeros(size(G)); use=referencePmp>0;
    normalization(use)=power(use)./referencePmp(use);
    currentScale=currentScale.*normalization;
else
    if isfield(inputs,'availablePower_W'), error('solarArrayIV:power','xtj_reference computes its own power; omit availablePower_W.'); end
end
label='NanoAvionics design proxy: XTJ curve shape, assumed wiring; NOT vendor verified';
result.time_s=t; result.irradiance_W_m2=G; result.cellTemperature_C=T;
result.Vmp_V=voltageScale*vm; result.Imp_A=currentScale.*jm;
result.Voc_V=voltageScale*voc; result.Isc_A=currentScale.*jsc;
result.Pmp_W=result.Vmp_V.*result.Imp_A;
result.testVoltage_V=p.testVoltage_V;
result.currentAtTestVoltage_A=[]; result.powerAtTestVoltage_W=[]; result.testVoltageAboveVoc=[];
if ~isempty(p.testVoltage_V)
    result.currentAtTestVoltage_A=currentScale.*currentAt(p.testVoltage_V/voltageScale+zeros(size(G)),voc,il,i0,a,rsRef);
    result.powerAtTestVoltage_W=p.testVoltage_V*result.currentAtTestVoltage_A;
    result.testVoltageAboveVoc=active & p.testVoltage_V>=result.Voc_V;
end
result.temperatureExtrapolated=active & (T<p.coefficientRange_C(1)|T>p.coefficientRange_C(2));
result.irradianceValidationUnavailable=active & abs(G-p.referenceIrradiance_W_m2)>1e-9;
result.referencePmp_W=referencePmp;
result.currentNormalization=normalization;
result.geometry=geometry; result.geometry.totalCells=p.totalCells;
result.geometry.cellArea_cm2=p.cellArea_cm2; result.geometry.cellArea_m2=p.cellArea_cm2/1e4;
result.wiringVerified=false;
result.verifiedForNanoAvionics=false; result.verifiedForFlightAssembly=false;
result.status=label; result.params=p;
ns=(1:p.totalCells).'; ns=ns(mod(p.totalCells,ns)==0); np=p.totalCells./ns;
result.topologyComparison=table(ns,np,ns*p.cellVmp_V,ns*p.cellVoc_V, ...
    np*p.cellArea_cm2*p.cellJmp_A_cm2, ...
    'VariableNames',{'SeriesCells','ParallelStrings','ReferenceVmp_V','ReferenceVoc_V','ReferenceImp_A'});
result.referenceFit=struct('diodeThermalVoltage_V',aRef,'seriesResistance_Ohm_cm2',rsRef, ...
    'lightCurrent_A_cm2',ilRef,'saturationCurrent_A_cm2',i0Ref,'maxRelativeResidual',fitResidual);
result.limitations=['Reduced single-diode fit, infinite shunt resistance; uniform illumination and temperature. ' ...
    'Wiring is an explicit scenario, not verified connector performance. In proxy mode current is normalized to the existing NanoAvionics power budget. ' ...
    'The panel boundary precedes external harness/PMAD losses. Bypass/shading/reverse bias are not simulated. ' ...
    'Temperature and irradiance response away from the calibration point need validation. ' ...
    'MPP and open circuit are distinct operating points; neither implies an actual connected voltage. Geometry footprint is treated as active area provisionally.'];
indices=find(active & result.Pmp_W>0);
if isempty(indices), selected=[]; else, selected=unique(indices(round(linspace(1,numel(indices),min(p.snapshotCount,numel(indices)))))); end
result.snapshotIndices=selected; result.curveVoltage_V=zeros(numel(selected),p.curvePoints); result.curveCurrent_A=result.curveVoltage_V;
result.snapshotSelectionNote=['Evenly spaced samples among generating time samples; ' ...
    'temperatures and irradiances are taken from the supplied orbit histories, not prescribed test conditions.'];
for k=1:numel(selected)
    j=selected(k); voltage=linspace(0,result.Voc_V(j),p.curvePoints);
    result.curveVoltage_V(k,:)=voltage;
    result.curveCurrent_A(k,:)=currentScale(j)*currentAt(voltage/voltageScale,voc(j),il(j),i0(j),a(j),rsRef);
end
if p.showSummary
    fprintf('\nSOLAR I-V: %s\n',label);
    fprintf('%gS%gP ASSUMED wiring; XTJ-only reference Vmp %.4f V, Voc %.4f V, Imp %.4f A, Pmp %.4f W.\n', ...
        p.seriesCells,p.parallelStrings,p.seriesCells*p.cellVmp_V,p.seriesCells*p.cellVoc_V, ...
        p.parallelStrings*p.cellArea_cm2*p.cellJmp_A_cm2, ...
        p.seriesCells*p.cellVmp_V*p.parallelStrings*p.cellArea_cm2*p.cellJmp_A_cm2);
    fprintf('Voltages follow specified series-cell count. Current scales with parallel strings.\n');
    if strcmp(p.mode,'power_consistent_proxy')
        fprintf('Proxy current is normalized to shared panel power; inspect currentNormalization and referencePmp_W.\n');
    end
    disp(result.topologyComparison);
    fprintf('Cell area %.6g cm^2; coverage %.4f%%; wiring and active area remain unverified.\n',p.cellArea_cm2,100*p.cellCoverage);
end
if p.makePlots
    fig=figure('WindowStyle','normal','Color','w','Name','Solar I-V estimate');
    finalPowerEPS.fitFigureToScreen(fig,p.figureSize_px);
    layout=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
    title(layout,sprintf('Solar I-V | %gS%gP assumed',p.seriesCells,p.parallelStrings));
    a=nexttile(layout); hold(a,'on'); b=nexttile(layout); hold(b,'on');
    colors=lines(max(1,numel(selected))); styles={'-','--',':','-.'};
    handles=gobjects(numel(selected),1);
    labels=cell(numel(selected),1);
    for k=1:numel(selected)
        j=selected(k); style=styles{1+mod(k-1,numel(styles))};
        handles(k)=plot(a,result.curveVoltage_V(k,:),result.curveCurrent_A(k,:), ...
            'Color',colors(k,:),'LineStyle',style,'LineWidth',p.plotLineWidth);
        plot(b,result.curveVoltage_V(k,:),result.curveVoltage_V(k,:).*result.curveCurrent_A(k,:), ...
            'Color',colors(k,:),'LineStyle',style,'LineWidth',p.plotLineWidth);
        labels{k}=sprintf('%.1f C, %.0f W/m^2',T(j),G(j));
    end
    grid(a,'on'); xlabel(a,'Array voltage [V]'); ylabel(a,'Current [A]'); title(a,'I-V curves');
    grid(b,'on'); xlabel(b,'Array voltage [V]'); ylabel(b,'Power [W]'); title(b,'P-V curves');
    if ~isempty(selected)
        voltageLimit=positiveLimit(result.curveVoltage_V,p.plotPaddingFraction/2);
        xlim(a,[0 voltageLimit]); xlim(b,[0 voltageLimit]);
    end
    ylim(a,[0 positiveLimit(result.curveCurrent_A,p.plotPaddingFraction)]);
    ylim(b,[0 positiveLimit(result.curveVoltage_V.*result.curveCurrent_A,p.plotPaddingFraction)]);
    if ~isempty(labels)
        lg=legend(a,handles,labels,'Orientation','horizontal','NumColumns',2,'Box','off');
        lg.Layout.Tile='south';
    end
    c=nexttile(layout); hold(c,'on');
    plot(c,t/60,result.Vmp_V,'LineWidth',p.plotLineWidth);
    plot(c,t/60,result.Voc_V,'--','LineWidth',p.plotLineWidth);
    grid(c,'on'); xlabel(c,'Time [min]'); ylabel(c,'Voltage [V]'); title(c,'Available voltages');
    legend(c,'MPP','Open circuit','Location','southoutside','Orientation','horizontal','Box','off');
    ylim(c,[0 positiveLimit([result.Vmp_V result.Voc_V],p.plotPaddingFraction)]);
    d=nexttile(layout); plot(d,t/60,result.Imp_A,'LineWidth',p.plotLineWidth); grid(d,'on');
    xlabel(d,'Time [min]'); ylabel(d,'MPP current [A]'); title(d,'Available current');
    ylim(d,[0 positiveLimit(result.Imp_A,p.plotPaddingFraction)]);
    if numel(t)>1, xlim(c,[t(1) t(end)]/60); xlim(d,[t(1) t(end)]/60); end
    set([a b c d],'FontSize',p.plotFontSize,'LineWidth',0.8);
end

end
function upper=positiveLimit(values,padding)
values=values(isfinite(values));
if isempty(values)||max(values)<=0, upper=1; else, upper=max(values)*(1+padding); end
end
function j=currentAt(v,voc,il,i0,a,rs)
lo=zeros(size(v+voc)); hi=voc+zeros(size(lo));
for k=1:60
    u=(lo+hi)/2; cur=il-i0.*expm1(u./a); value=u-cur*rs;
    up=value<v; lo(up)=u(up); hi(~up)=u(~up);
end
j=max(0,il-i0.*expm1((lo+hi)/2./a)); j(v>=voc)=0;
end
function x=signal(x,n)
validateattributes(x,{'numeric'},{'real','finite','vector','nonempty'}); x=x(:).';
if isscalar(x), x=repmat(x,1,n); end
if numel(x)~=n, error('solarArrayIV:length','Signals must be scalar or match time.'); end
end
function [a,rs,il,i0,err]=fitReference(p)
% Fit Voc, Jsc, Vmp, Jmp and d(V*J)/dV=0 at MPP, with Rsh=infinity.
persistent cachedKey cachedFit
key=[p.cellVoc_V,p.cellVmp_V,p.cellJsc_A_cm2,p.cellJmp_A_cm2];
if isequal(key,cachedKey)
    a=cachedFit(1); rs=cachedFit(2); il=cachedFit(3); i0=cachedFit(4); err=cachedFit(5); return
end
options=optimset('Display','off','TolX',1e-11,'TolFun',1e-18,'MaxIter',2000,'MaxFunEvals',5000);
z=fminsearch(@objective,log([0.067 4.3]),options); [r,a,rs,il,i0]=residual(z); err=max(abs(r));
if err>1e-7, error('solarArrayIV:fit','Reference points do not admit the configured reduced diode model.'); end
cachedKey=key; cachedFit=[a,rs,il,i0,err];
    function f=objective(z)
        r=residual(z); f=sum(r.^2);
    end
    function [r,a,rs,il,i0]=residual(z)
        a=exp(z(1)); rs=exp(z(2));
        v=p.cellVoc_V; sc=p.cellJsc_A_cm2; m=p.cellJmp_A_cm2; vm=p.cellVmp_V;
        il=NaN; i0=NaN;
        if ~isfinite(a+rs)||a<0.005||a>1||sc*rs>=v||v/a>650
            r=[1e3;1e3]; return;
        end
        i0=sc/(exp(v/a)-exp(sc*rs/a)); il=i0*expm1(v/a);
        predicted=il-i0*expm1((vm+m*rs)/a);
        r=[(predicted-m)/sc;(a/(il-m+i0)+rs-vm/m)/(vm/m)];
        if any(~isfinite(r)), r=[1e3;1e3]; end
    end
end
