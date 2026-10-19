%-------------------------------------------------------------------------%
%                                                                         %
%  testj77.m                                                              %
%                                                                         %
%  Written by:  David L. Huestis, Molecular Physics Laboratory            %
%                                                                         %
%  Copyright (c) 1999  SRI International                                  %
%  All Rights Reserved                                                    %
%                                                                         %
%  This software is provided on an as is basis; without any               %
%  warranty; without the implied warranty of merchantability or           %
%  fitness for a particular purpose.                                      %
%                                                                         %
% Last modified:   2015/08/26   M. Mahooti                                %
%                                                                         %
%-------------------------------------------------------------------------%
clc
clear all
format long g

mz = 2500;
% Tinf = 600.0;
Tinf = 1000.0;
z = (0:mz);
v = zeros(7,1);
[T, CN2, CO2, CO, CAr, CHe, CH, CM, WM] = j77sri(mz, Tinf, z);

fileID = fopen('t1000.out','w');

for i = 0:mz
    if( i <= 80 )
        if( mod(i,5) ~= 0 )
            continue;
        end
    elseif( i <= 100 )
            v(1) = CN2(i+1);
            v(2) = CO2(i+1);
            v(3) = CO(i+1);
            v(4) = CAr(i+1);
            v(5) = CHe(i+1);
            v(6) = CH(i+1);
            v(7) = CM(i+1);
            for j=1:7
                if( v(j) > 1.26e-16 )
                    v(j) = log10( v(j) ) + 6.0;
                else
                    v(j) = -9.9;
                end
            end
            fprintf(fileID,'%5i %8.2f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %7.3f\n',...
                i,T(i+1),v(1),v(2),v(3),v(4),v(5),v(6),v(7),WM(i+1));
            continue;
    elseif( i <= 110 )
        if( mod(i,2) ~= 0 )
            continue;
        end
    elseif( i <= 160 )
        if( mod(i,5) ~= 0 )
            continue;
        end
    elseif( i <= 400 )
        if( mod(i,10) ~= 0 )
            continue;
        end
    elseif( i <= 1000 )
        if( mod(i,20) ~= 0 )
            continue;
        end
    elseif( i <= 1500 )
        if( mod(i,50) ~= 0 )
            continue;
        end
    elseif( mod(i,100) ~= 0 )
        continue;
    end
    v(1) = CN2(i+1);
    v(2) = CO2(i+1);
    v(3) = CO(i+1);
    v(4) = CAr(i+1);
    v(5) = CHe(i+1);
    v(6) = CH(i+1);
    v(7) = CM(i+1);
    for j=1:7
        if( v(j) > 1.26e-16 )
            v(j) = log10( v(j) ) + 6.0;
        else
            v(j) = -9.9;
        end
    end
    fprintf(fileID,'%5i %8.2f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %7.3f\n',...
        i,T(i+1),v(1),v(2),v(3),v(4),v(5),v(6),v(7),WM(i+1));
end
fclose(fileID);


