function [Ftgt] = matchedfilter_psf(npix)
%Generate Event Match Filter PSF
% %% Telescope Parameters
%     %%Officina Stellare RH200
%     D=0.2;%m telescope diameter
%     f=0.6;%m focal length
%     N=f/D;
%% Lens Parameters
    D=0.1;%m telescope diameter
    f=0.475;%m focal length
    N=f/D;
%% Detector Parameters
%Prophesee EVK3 Gen4
l1=400e-9;
l2=1000e-9;
l=(l1+l2)/2;
p=4.86e-6; %(m) Pixel Size

% npix=50;

ifov=2*atan(p/(2*N*D));


platescale=1/f;%arcsec/m
% return


%PSF variance
sigsq=10*(0.45*l*N)^2*platescale^2; %Replace with sigma due to seeing (6-7 arcsec in Boulder per Rose Smith)

% return
%% Calculate Fraction of Total Flux in Each Pixel

mux=npix*ifov/2;
muy=npix*ifov/2;

x= (0:npix)*ifov;
y= (0:npix)*ifov;
[X,Y]=meshgrid(x,y);

Xm=(X(1:end-1,1:end-1)-mux)/sqrt(2*sigsq);
Xp=(X(1:end-1,2:end)-mux)/sqrt(2*sigsq);
Ym=(Y(1:end-1,1:end-1)-muy)/sqrt(2*sigsq);
Yp=(Y(2:end,1:end-1)-muy)/sqrt(2*sigsq);
Ftgt=1/4*(erf(Xp)-erf(Xm)).*(erf(Yp)-erf(Ym)); % See OneNote Derivation
end

