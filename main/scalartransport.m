%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% scalartransport.m
% usage: calculate scalar transport in swe2d.m
% using First-order upwinding and theta method 
% for velocity( no diffusion here!)
%
% Yun Zhang @ Stanford
% 4/4/2013 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function C_new = scalartransport(u_old,v_old,u_new,v_new,xc,yc,xu,yu,xv,yv,C,Huwx,Huwy,...
                                 dx,dy,dt,n,theta,h,d,H_old,L,W,rhs_C,boundaryregion_type2,...
                                 Cb_type2,boundaryregion_type3,Cb_type3)

global cellmark 

getvariables

Cuwx=zeros(size(u_new));
Cuwy=zeros(size(v_new));

H_new = h+d;

utheta=theta*u_new+(1-theta)*u_old;
vtheta=theta*v_new+(1-theta)*v_old;

% First-order upwind values on the u faces
Cuwx(2:end-1,:)=(utheta(2:end-1,:)>0).*C(1:end-1,:)+(utheta(2:end-1,:)<=0).*C(2:end,:);

% First-order upwind values on the v faces
Cuwy(:,2:end-1)=(vtheta(:,2:end-1)>0).*C(:,1:end-1)+(vtheta(:,2:end-1)<=0).*C(:,2:end);


% This sets the value at the face regardless of whether upwind or
% downwind.  Need to fix so that it updates based on upwind,
% i.e. uses interior value instead of boundary value...
[Cuwx,Cuwy]=boundaryconditions_face(C,u_old,v_old,Cuwx,Cuwy,Huwx*dy,Huwy*dx,xu,yu,xv,yv,n,boundaryregion_type2,[],[],[],Cb_type2,'C');
C=boundaryconditions_cell(C,xc,yc,n,boundaryregion_type3,Cb_type3);

if(isempty(rhs_C))
    rhs_C=zeros(size(C));
end

C_new=zeros(size(C));
C_new=1./H_new.*(H_old.*C + dt*rhs_C ...
                 -dt*(Huwx(2:end,:).*utheta(2:end,:).*Cuwx(2:end,:) ...
                      -Huwx(1:end-1,:).*utheta(1:end-1,:).*Cuwx(1:end-1,:))/dx ...
                 -dt*(Huwy(:,2:end).*vtheta(:,2:end).*Cuwy(:,2:end) ...
                      -Huwy(:,1:end-1).*vtheta(:,1:end-1).*Cuwy(:,1:end-1))/dy);

% Don't touch the type 3 boundary cells
C_new(cellmark==3)=C(cellmark==3);
C_new(H_new==0)=0;

% Need to set inactive cells to zero since they may have been
% updated at boundary edges that are inside the domain.
C_new(cellmark==0)=0;

% Compute fluxes here for assignments 2 and 3
if(exist('FLUXES','var') & FLUXES)
    fluxes;
end
