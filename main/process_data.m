% 
% This script is called at the end of each time step. It can be
% used to access/process/output data for analysis.  Variables
% include:
%
% xc : Cell-centered x locations in m (size Ni x Nj)
% yc : Cell-centered y locations in m (size Ni x Nj)
% u : Face-centered u velocity in m/s (size Ni+1 x Nj)
% v : Face-centered v velocity in m/s (size Ni x Nj+1)
% C : Cell-centered SSC or scalar in mg/l (size Ni x Nj)
% 

in=find(cellmark~=0);
H=(h+d);

% This is the volume in the domain at time step n
Volume(n)=sum(sum(H(in)))*dx*dy;
  
% This is the maximum concentration at time step n
Cmax(n)=max(max(C(cellmark~=0)));

% This is the grid cell coordinate where the maximum concentration
% exists; Use the first index (imax(1),jmax(1)) in case there are 
% repeated values
[imax,jmax]=find(C==Cmax(n));
xc_imax = xc(imax(1),jmax(1)); 
yc_imax = yc(imax(1),jmax(1));

%fprintf('Cmax is %.2f mg/l located at (%.2f, %.2f).\n',...
%        Cmax(n),xc_imax,yc_imax);