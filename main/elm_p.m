%
% Compute the velocity vectors at the location of the Lagrangian
% tracebacks.
%
% Oliver Fringer
% Yun Zhang
% Stanford University
%
function [uminus,vminus]=elm_p(xu,yu,xv,yv,xc,yc,d,H_small,u,v,u_vface,v_uface,dt,xpoly,ypoly)

global edge_u_ij edge_u_ijm1 poly_ij edge_v_ij edge_v_im1j markuface ...
    markvface

dx = xu(2,1)-xu(1,1);
dy = yu(1,2)-yu(1,1);
[Ni,Nj]=size(d);

[xp,yp]=ndgrid([0:Ni]*dx,[0:Nj]*dy);
up = zeros(Ni+1,Nj+1);
vp = zeros(Ni+1,Nj+1);
up(:,2:Nj) = 0.5*(u(:,1:Nj-1)+u(:,2:Nj));
vp(2:Ni,:) = 0.5*(v(1:Ni-1,:)+v(2:Ni,:));

% Use the nonzero neighboring value of the velocity to estimate the
% velocity on the corners.
if 1
    a=u(edge_u_ij);
    b=u(edge_u_ijm1);
    up(poly_ij)=0;(a==b).*a+(a~=b).*((a==0).*b+(b==0).*a);

    a=v(edge_v_ij);
    b=v(edge_v_im1j);
    vp(poly_ij)=0;(a==b).*a+(a~=b).*((a==0).*b+(b==0).*a);
end

u_vface = 0.5*(up(1:Ni,:)+up(2:Ni+1,:));
v_uface = 0.5*(vp(:,1:Nj)+up(:,2:Nj+1));

% Calculate the Lagrangian traceback from the u-faces
xminus_uface = xu - dt*u;
yminus_uface = yu - dt*v_uface;

% Only update the velocities for which the traceback is
% within the domain, or where the depth is greater than H_small
%in = inpolygon(xminus_uface,yminus_uface,xpoly,ypoly);
uminus=u;
uminus(markuface==0) = interp2(xp',yp',up',...
                               xminus_uface(markuface==0),...
                               yminus_uface(markuface==0),'linear');
uminus(isnan(uminus))=u(isnan(uminus));

% Calculate the Lagrangian traceback from the v-faces
xminus_vface = xv - dt*u_vface;
yminus_vface = yv - dt*v;

% Only update the velocities for which the traceback is
% within the domain, or where the depth is greater than H_small
%in = inpolygon(xminus_vface,yminus_vface,xpoly,ypoly);
vminus=v;
vminus(markvface==0) = interp2(xp',yp',vp',...
                               xminus_vface(markvface==0),...
                               yminus_vface(markvface==0),'linear');
vminus(isnan(vminus))=v(isnan(vminus));

% Don't do ELM on boundary edges.
uminus(markuface>0)=u(markuface>0);
vminus(markvface>0)=v(markvface>0);

