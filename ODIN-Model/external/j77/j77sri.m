%-------------------------------------------------------------------------%
%                                                                         %
%  j77sri.m                                                               %
%                                                                         %
%  Written by:  David L. Huestis, Molecular Physics Laboratory            %
%                                                                         %
%  Copyright (c) 1999,2002  SRI International                             %
%  All Rights Reserved                                                    %
%                                                                         %
%  This software is provided on an as is basis; without any               %
%  warranty; without the implied warranty of merchantability or           %
%  fitness for a particular purpose.                                      %
%                                                                         %
%-------------------------------------------------------------------------%
%                                                                         %
%      Given an exospheric temperature, this subroutine returns model     %
%      atmospheric altitude profiles of temperature, the number           %
%      densities of N2, O2, O, Ar, He, H, the sum thereof, and the        %
%      molecular weight.                                                  %
%                                                                         %
%      For altitudes of 90 km and above, we use the 1977 model of         %
%      Jacchia [Ja77].  H-atom densities are returned as non-zero         %
%      for altitudes of 150 km and above if the maximum altitude          %
%      requested is 500 km or more.                                       %
%                                                                         %
%      For altitudes of 85 km and below we use the 1976 U. S. Standard    %
%      Atmosphere, as coded by Carmichael [Ca99], which agrees with       %
%      Table III.1 (pp 422-423) of Chamberlain and Hunten [CH87]          %
%	   and Table I (pp 50-73) of the "official" U.S. Standard             %
%	   Atmosphere 1976 [COESA76].                                         %
%                                                                         %
%      For altitudes from 86 to 89 km we calculate the extent of          %
%      oxygen dissociation and the effective molecular weight by a        %
%      polynomial fit connecting the O-atom mole fraction at 86 km        %
%      from Chamberlain and Hunten Table III.4 (p 425) [CH87] and         %
%      the O-atmom mole fractions at 90, 91, and 92 km from Jacchia       %
%      1977 [Ja77] for an exospheric temperature of 1000 K.  For          %
%      graphical continunity, the same formulas are used to calculate     %
%      O-atom densities for altitudes of 85 km and below.                 %
%                                                                         %
%  USAGE:                                                                 %
%              program main                                               %
%              integer maxz    ! INPUT:  highest altitude (km)            %
%              parameter (maxz=2500)   ! for example                      %
%              real Tinf,      ! INPUT:  exospheric temp (K)              %
%           *    Z(0:maxz),    ! OUTPUT: altitude (km)                    %
%           *    T(0:maxz),    ! OUTPUT: temperature (K)                  %
%           *    CN2(0:maxz),  ! OUTPUT: [N2] (1/cc)                      %
%           *    CO2(0:maxz),  ! OUTPUT: [O2] (1/cc)                      %
%           *    CO(0:maxz),   ! OUTPUT: [O] (1/cc)                       %
%           *    CAr(0:maxz),  ! OUTPUT: [Ar] (1/cc)                      %
%           *    CHe(0:maxz),  ! OUTPUT: [He] (1/cc)                      %
%           *    CH(0:maxz),   ! OUTPUT: [H] (1/cc)                       %
%           *    CM(0:maxz),   ! OUTPUT: [M] (1/cc)                       %
%           *    WM(0:maxz)    ! OUTPUT: molecular weight (g)             %
%              call j77sri(maxz,Tinf,Z,T,CN2,CO2,CO,CAr,CHe,CH,CM,WM)     %
%              end                                                        %
%                                                                         %
%  REFERENCES:                                                            %
%                                                                         %
%      Ca99    R. Carmichael, "Fortran (90) coding of Atmosphere,"        %
%              (http://www.pdas.com/atmosf90.htm, March 1, 1999).         %
%                                                                         %
%      CH87    J. W. Chamberlain and D. M. Hunten, "Theory of             %
%              Planetary Atmospheres," (Academic Press, NY, 1987).        %
%                                                                         %
%	COESA76	U.S. Committee on Extension to the Standard                   %
%		Atmosphere, "U.S. Standard Atmospheres 1976"                      %
%		(USGPO, Washington, DC, 1976).                                    %
%                                                                         %
%      Ja77    L. G. Jacchia, "Thermospheric Temperature, Density         %
%              and Composition: New Models," SAO Special Report No.       %
%              375 (Smithsonian Institution Astrophysical                 %
%              Observatory, Cambridge, MA, March 15, 1977).               %
%                                                                         %
% Last modified:   2015/08/26   M. Mahooti                                %
%                                                                         %
%-------------------------------------------------------------------------%
function [T, CN2, CO2, CO, CAr, CHe, CH, CM, WM] = j77sri(maxz, Tinf, Z)

pi2 = 1.57079632679;
wm0 = 28.96; wmN2 = 28.0134; wmO2 = 31.9988; wmO = 15.9994; wmAr = 39.948;
wmHe = 4.0026; wmH = 1.0079;
qN2 = 0.78110; qO2 = 0.20955; qAr = 0.009343; qHe = 0.000005242;
T = (0:maxz); CN2 = (0:maxz); CO2 = (0:maxz); CO = (0:maxz);
CAr = (0:maxz); CHe = (0:maxz); CH = (0:maxz); CM = (0:maxz);
WM = (0:maxz); E5M = (0:10); E6P = (0:10);

for iz = 0:maxz
    Z(iz+1) = iz;
    CH(iz+1) = 0;
    if( iz <= 85 )
        h = Z(iz+1)*6369.0/(Z(iz+1)+6369.0);
        if( iz <= 32 )
            if( iz <= 11 )
            hbase = 0.0;
            pbase = 1.0;
            tbase = 288.15;
            tgrad = -6.5;
            T(iz+1) = tbase + tgrad*(h-hbase);
            x = (tbase/T(iz+1))^(34.163195/tgrad);
            CM(iz+1) = 2.547e19*(288.15/T(iz+1))*pbase*x;
            y = 10.0^(-3.7469+(iz-85)*(0.226434-(iz-85)*5.945e-3));
            x = 1 - y;
            WM(iz+1) = wm0*x;
            CN2(iz+1) = qN2*CM(iz+1);
            CO(iz+1) = 2.0*y*CM(iz+1);
            CO2(iz+1) = (x*qO2-y)*CM(iz+1);
            CAr(iz+1) = qAr*CM(iz+1);
            CHe(iz+1) = qHe*CM(iz+1);
            CH(iz+1) = 0;
            continue;
            elseif( iz <= 20 )
                hbase = 11;
                pbase = 2.233611e-1;
                tbase = 216.65;
                tgrad = 0;
                T(iz+1) = tbase;
                x = exp(-34.163195*(h-hbase)/tbase);
                CM(iz+1) = 2.547e19*(288.15/T(iz+1))*pbase*x;
                y = 10.0^(-3.7469+(iz-85)*(0.226434-(iz-85)*5.945e-3));
                x = 1 - y;
                WM(iz+1) = wm0*x;
                CN2(iz+1) = qN2*CM(iz+1);
                CO(iz+1) = 2.0*y*CM(iz+1);
                CO2(iz+1) = (x*qO2-y)*CM(iz+1);
                CAr(iz+1) = qAr*CM(iz+1);
                CHe(iz+1) = qHe*CM(iz+1);
                CH(iz+1) = 0;
                continue;
            else
                hbase = 20.0;
                pbase = 5.403295e-2;
                tbase = 216.65;
                tgrad = 1;
                T(iz+1) = tbase + tgrad*(h-hbase);
                x = (tbase/T(iz+1))^(34.163195/tgrad);
                CM(iz+1) = 2.547e19*(288.15/T(iz+1))*pbase*x;
                y = 10.0^(-3.7469+(iz-85)*(0.226434-(iz-85)*5.945e-3));
                x = 1 - y;
                WM(iz+1) = wm0*x;
                CN2(iz+1) = qN2*CM(iz+1);
                CO(iz+1) = 2.0*y*CM(iz+1);
                CO2(iz+1) = (x*qO2-y)*CM(iz+1);
                CAr(iz+1) = qAr*CM(iz+1);
                CHe(iz+1) = qHe*CM(iz+1);
                CH(iz+1) = 0;
                continue;
            end
        elseif( iz <= 51 )
            if( iz <= 47 )
                hbase = 32.0;
                pbase = 8.5666784e-3;
                tbase = 228.65;
                tgrad = 2.8;
                T(iz+1) = tbase + tgrad*(h-hbase);
                x = (tbase/T(iz+1))^(34.163195/tgrad);
                CM(iz+1) = 2.547e19*(288.15/T(iz+1))*pbase*x;
                y = 10.0^(-3.7469+(iz-85)*(0.226434-(iz-85)*5.945e-3));
                x = 1 - y;
                WM(iz+1) = wm0*x;
                CN2(iz+1) = qN2*CM(iz+1);
                CO(iz+1) = 2.0*y*CM(iz+1);
                CO2(iz+1) = (x*qO2-y)*CM(iz+1);
                CAr(iz+1) = qAr*CM(iz+1);
                CHe(iz+1) = qHe*CM(iz+1);
                CH(iz+1) = 0;
                continue;
            else
                hbase = 47;
                pbase = 1.0945601e-3;
                tbase = 270.65;
                tgrad = 0;
                T(iz+1) = tbase;
                x = exp(-34.163195*(h-hbase)/tbase);
                CM(iz+1) = 2.547e19*(288.15/T(iz+1))*pbase*x;
                y = 10.0^(-3.7469+(iz-85)*(0.226434-(iz-85)*5.945e-3));
                x = 1 - y;
                WM(iz+1) = wm0*x;
                CN2(iz+1) = qN2*CM(iz+1);
                CO(iz+1) = 2.0*y*CM(iz+1);
                CO2(iz+1) = (x*qO2-y)*CM(iz+1);
                CAr(iz+1) = qAr*CM(iz+1);
                CHe(iz+1) = qHe*CM(iz+1);
                CH(iz+1) = 0;
                continue;
            end
        elseif( iz <= 71 )
            hbase = 51.0;
            pbase = 6.6063531e-4;
            tbase = 270.65;
            tgrad = -2.8;
            T(iz+1) = tbase + tgrad*(h-hbase);
            x = (tbase/T(iz+1))^(34.163195/tgrad);
            CM(iz+1) = 2.547e19*(288.15/T(iz+1))*pbase*x;
            y = 10.0^(-3.7469+(iz-85)*(0.226434-(iz-85)*5.945e-3));
            x = 1 - y;
            WM(iz+1) = wm0*x;
            CN2(iz+1) = qN2*CM(iz+1);
            CO(iz+1) = 2.0*y*CM(iz+1);
            CO2(iz+1) = (x*qO2-y)*CM(iz+1);
            CAr(iz+1) = qAr*CM(iz+1);
            CHe(iz+1) = qHe*CM(iz+1);
            CH(iz+1) = 0;
            continue;
        else
            hbase = 71.0;
            pbase = 3.9046834e-5;
            tbase = 214.65;
            tgrad = -2.0;
            T(iz+1) = tbase + tgrad*(h-hbase);
            x = (tbase/T(iz+1))^(34.163195/tgrad);
            CM(iz+1) = 2.547e19*(288.15/T(iz+1))*pbase*x;
            y = 10.0^(-3.7469+(iz-85)*(0.226434-(iz-85)*5.945e-3));
            x = 1 - y;
            WM(iz+1) = wm0*x;
            CN2(iz+1) = qN2*CM(iz+1);
            CO(iz+1) = 2.0*y*CM(iz+1);
            CO2(iz+1) = (x*qO2-y)*CM(iz+1);
            CAr(iz+1) = qAr*CM(iz+1);
            CHe(iz+1) = qHe*CM(iz+1);
            CH(iz+1) = 0;
            continue;
        end
    elseif( iz <= 89 )
        T(iz+1) = 188.0;
        y = 10.0^(-3.7469+(iz-85)*(0.226434-(iz-85)*5.945e-3));
        WM(iz+1) = wm0*(1-y);
        CM(iz+1) = CM(iz)*(T(iz)/T(iz+1))*(WM(iz+1)/WM(iz))...
            * exp( - 0.5897446*( ...
            (WM(iz)/T(iz))*(1+Z(iz)/6356.766)^(-2)...
            + (WM(iz+1)/T(iz+1))*(1+Z(iz+1)/6356.766)^(-2) ));
        y = 10.0^(-3.7469+(iz-85)*(0.226434-(iz-85)*5.945e-3));
        x = 1 - y;
        WM(iz+1) = wm0*x;
        CN2(iz+1) = qN2*CM(iz+1);
        CO(iz+1) = 2.0*y*CM(iz+1);
        CO2(iz+1) = (x*qO2-y)*CM(iz+1);
        CAr(iz+1) = qAr*CM(iz+1);
        CHe(iz+1) = qHe*CM(iz+1);
        CH(iz+1) = 0;
        continue;
    elseif( iz <= 90 )
        T(iz+1) = 188.0;
    elseif( Tinf < 188.1 )
        T(iz+1) = 188.0;
    else
        x = 0.0045 * (Tinf-188.0);
        Tx = 188.0 + 110.5 * log( x + sqrt(x*x+1.0) );
        Gx = pi2*1.9*(Tx - 188.0)/(125.0-90.0);
        if( iz <= 125 )
            T(iz+1) = Tx + ((Tx-188.0)/pi2)...
                * atan( (Gx/(Tx-188.0))*(Z(iz+1)-125.0)...
                * (1.0 + 1.7*((Z(iz+1)-125.0)/(Z(iz+1)-90.0))^2));
        else
            T(iz+1) = Tx + ((Tinf-Tx)/pi2)...
                * atan( (Gx/(Tinf-Tx))*(Z(iz+1)-125.0)...
                * (1.0 + 5.5e-5*(Z(iz+1)-125.0)^2));
        end
    end
    if( iz <= 100 )
        x = iz - 90;
        E5M(iz+1-90) = 28.89122 + x*(-2.83071e-2 ...
            + x*(-6.59924e-3 + x*(-3.39574e-4 ...
            + x*(+6.19256e-5 + x*(-1.84796e-6) ))));
        if( iz <= 90 )
            E6P(1) = 7.145e13*T(91);
        else
            G0 = (1+Z(iz)/6356.766)^(-2);
            G1 = (1+Z(iz+1)/6356.766)^(-2);
            E6P(iz+1-90) = E6P(iz+1-91)*exp( - 0.5897446*( ...
                G1*E5M(iz+1-90)/T(iz+1) + G0*E5M(iz+1-91)/T(iz) ) );
        end
        x = E5M(iz+1-90)/wm0;
        y = E6P(iz+1-90)/T(iz+1);
        CN2(iz+1) = qN2*y*x;
        CO(iz+1) = 2.0*(1.0 - x)*y;
        CO2(iz+1) = (x*(1.0+qO2)-1.0)*y;
        CAr(iz+1) = qAr*y*x;
        CHe(iz+1) = qHe*y*x;
        CH(iz+1) = 0;
    else
        G0 = (1+Z(iz)/6356.766)^(-2);
        G1 = (1+Z(iz+1)/6356.766)^(-2);
        x =  0.5897446*( G1/T(iz+1) + G0/T(iz) );
        y = T(iz)/T(iz+1);
        CN2(iz+1) = CN2(iz)*y*exp(-wmN2*x);
        CO2(iz+1) = CO2(iz)*y*exp(-wmO2*x);
        CO(iz+1)  = CO(iz)*y*exp(-wmO*x);
        CAr(iz+1) = CAr(iz)*y*exp(-wmAr*x);
        CHe(iz+1) = CHe(iz)*(y^0.62)*exp(-wmHe*x);
        CH(iz+1) = 0;
    end
    continue;
%     y = 10.0^(-3.7469+(iz+1-85)*(0.226434-(iz+1-85)*5.945E-3));
%           x = 1 - y;
%           WM(iz+1) = wm0*x
%           CN2(iz+1) = qN2*CM(iz+1)
%           CO(iz+1) = 2.0*y*CM(iz+1)
%           CO2(iz+1) = (x*qO2-y)*CM(iz+1)
%           CAr(iz+1) = qAr*CM(iz+1)
%           CHe(iz+1) = qHe*CM(iz+1)
%           CH(iz+1) = 0
end
%  ----------------------------------------------------------------------
%       Add Jacchia 1977 empirical corrections to [O] and [O2]
%  ----------------------------------------------------------------------
for iz = 90:maxz
    CO2(iz+1) = CO2(iz+1)...
        *( 10.0^(-0.07*(1.0+tanh(0.18*(Z(iz+1)-111.0)))) );
    CO(iz+1) = CO(iz+1)...
        *( 10.0^(-0.24*exp(-0.009*(Z(iz+1)-97.7)^2)) );
    CM(iz+1) = CN2(iz+1)+CO2(iz+1)+CO(iz+1)+CAr(iz+1)+CHe(iz+1)+CH(iz+1);
    WM(iz+1) = ( wmN2*CN2(iz+1)+wmO2*CO2(iz+1)+wmO*CO(iz+1)...
        +wmAr*CAr(iz+1)+wmHe*CHe(iz+1)+wmH*CH(iz+1) ) / CM(iz+1);
end
%  ----------------------------------------------------------------------
%       Calculate [H] from Jacchia 1997 formulas if maxz >= 500.
%  ----------------------------------------------------------------------
if( maxz >= 500 )
    phid00 = 10.0^( 6.9 + 28.9*Tinf^(-0.25) ) / 2.e20;
    phid00 = phid00 * 5.24e2;
    H_500 = 10.0^( -0.06 + 28.9*Tinf^(-0.25) );
    for iz = 150:maxz
        phid0 = phid00/sqrt(T(iz+1));
        WM(iz+1) = wmH*0.5897446*( (1.0+Z(iz+1)/6356.766)^(-2) )...
            / T(iz+1) + phid0;
        CM(iz+1) = CM(iz+1)*phid0;
    end
    y = WM(151);
    WM(151) = 0;
    for iz = 151:maxz
        x = WM(iz) + (y+WM(iz+1));
        y = WM(iz+1);
        WM(iz+1) = x;
    end
    for iz = 150:maxz
        WM(iz+1) = exp( WM(iz+1) ) * ( T(iz+1)/T(151) )^0.75;
        CM(iz+1) = WM(iz+1)*CM(iz+1);
    end
    y = CM(151);
    CM(151) = 0;
    for iz = 151:maxz
        x = CM(iz) + 0.5*(y+CM(iz+1));
        y = CM(iz+1);
        CM(iz+1) = x;
    end
    for iz = 150:maxz
        CH(iz+1) = ( WM(501)/WM(iz+1) ) * (H_500 - (CM(iz+1)-CM(501)) );
    end
    for iz=150:maxz
        CM(iz+1) = CN2(iz+1)+CO2(iz+1)+CO(iz+1)+CAr(iz+1)+CHe(iz+1)+CH(iz+1);
        WM(iz+1) = ( wmN2*CN2(iz+1)+wmO2*CO2(iz+1)+wmO*CO(iz+1)...
            +wmAr*CAr(iz+1)+wmHe*CHe(iz+1)+wmH*CH(iz+1) ) / CM(iz+1);
    end
end

