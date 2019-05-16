function zbnew = exner(u,v,zb,dx,dy,dt,H,tan_deltas,C,rhs_C,n)

getvariables;

% Grid size
[Ni,Nj] = size(zb);

% Need to return zbnew given zb and rhs_C
zbnew = zeros(Ni,Nj);

