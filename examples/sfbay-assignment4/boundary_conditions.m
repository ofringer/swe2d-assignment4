Tday=86400;
t=[1:nsteps]*dt; 

% Simplified tide at the ocean boundary that ramps up over one day
% with the function 1 - exp(-t/Tday)
amp=0.5;
TM2=12*3600; % Use 12 hours to ensure periodicity with Tday
omegaM2=2*pi/TM2;
hb=hoffset+TIDES*amp*sin(omegaM2*t).*(1-exp(-t/Tday)); 

% This polygon surrounds the open ocean boundary
boundaryregion_type3{1} = 1e4*...
    [0.1085    4.0875
     0.9669    4.1366
     0.9423    4.4309
     0.1085    4.4677];
Cb_type3{1} = 0;
hb_type3{1} = hb;

boundaryregion_type2 = [];
Qb_type2 = [];
Cb_type2 = [];

