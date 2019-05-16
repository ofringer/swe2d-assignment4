function [Cdx,Cdy]=dragpatch(Cdx,Cdy,xu,yu,xv,yv)

xbox = [600 700 700 600 600];
ybox = [200 200 300 300 200];

inu=inpolygon(xu,yu,xbox,ybox);
inv=inpolygon(xv,yv,xbox,ybox);

Cdveg = 1;
Cdx(inu) = Cdveg;
Cdy(inv) = Cdveg;
