%
% Determine whether a set of points lie within the domain by
% returning true for points in which the depth is nonzero.
%
function flag = incheck(xm,ym,x,y,d,Hmin)

dm = interp2(x,y,d,xm,ym);
flag = true(size(xm));
flag(find(dm<Hmin | isnan(dm)))=false;

