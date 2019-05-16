function [u,v,Cu,Cv] = advection_ghost(u,v,dx,dy,Cu,Cv,n,dt,nu)

global_pointers;

[Ni,Nj]=size(cellmark);

ug = zeros(Ni+3,Nj+2);
vg = zeros(Ni+2,Nj+3);

ug(2:end-1,2:end-1) = u;
vg(2:end-1,2:end-1) = v;

[is,js]=find(cellmark~=0);
for m=1:length(is)
    if(is(m)==1)
        u(is(m),js(m))=u(is(m)+1,js(m));
    elseif(is(m)==Ni)
        u(is(m)+1,js(m))=u(is(m),js(m));
    elseif(cellmark(is(m),js(m))~=cellmark(is(m)-1,js(m)))
        u(is(m),js(m))=u(is(m)+1,js(m));        
    else
        u(is(m)+1,js(m))=u(is(m),js(m));        
    end        

    if(js(m)==1)
        v(is(m),js(m))=v(is(m),js(m)+1);
    elseif(js(m)==Nj)
        v(is(m),js(m)+1)=v(is(m),js(m));
    elseif(cellmark(is(m),js(m))~=cellmark(is(m),js(m)-1))
        v(is(m),js(m))=v(is(m),js(m)+1);        
    else
        v(is(m),js(m)+1)=v(is(m),js(m));        
    end        
end

up = zeros(Ni+1,Nj+1);
vp = zeros(Ni+1,Nj+1);
up(:,2:Nj) = 0.5*(u(:,1:Nj-1)+u(:,2:Nj));
vp(2:Ni,:) = 0.5*(v(1:Ni-1,:)+v(2:Ni,:));
        
% Use the nonzero neighboring value of the velocity to estimate the
% velocity on the corners.
a=u(edge_u_ij);
b=u(edge_u_ijm1);
% Free slip
%up(poly_ij)=(a==b).*a+(a~=b).*((a==0).*b+(b==0).*a);

a=v(edge_v_ij);
b=v(edge_v_im1j);
% Free slip
%vp(poly_ij)=(a==b).*a+(a~=b).*((a==0).*b+(b==0).*a);

ux = zeros(Ni+1,Nj+1);
uy = zeros(Ni+1,Nj+1);
vx = zeros(Ni+1,Nj+1);
vy = zeros(Ni+1,Nj+1);

uxx = zeros(Ni+1,Nj);
uyy = zeros(Ni+1,Nj);
vxx = zeros(Ni,Nj+1);
vyy = zeros(Ni,Nj+1);

ux(2:end-1,2:end-1) = (up(3:end,2:end-1)-up(1:end-2,2:end-1))/(2*dx);
uy(2:end-1,2:end-1) = (up(2:end-1,3:end)-up(2:end-1,1:end-2))/(2*dy);
uxx(2:end-1,:) = (u(3:end,:)-2*u(2:end-1,:)+u(1:end-2,:))/dx^2;
uyy(:,2:end-1) = (u(:,3:end)-2*u(:,2:end-1)+u(:,1:end-2))/dy^2;

vx(2:end-1,2:end-1) = (vp(3:end,2:end-1)-vp(1:end-2,2:end-1))/(2*dx);
vy(2:end-1,2:end-1) = (vp(2:end-1,3:end)-vp(2:end-1,1:end-2))/(2*dy);
vxx(2:end-1,:) = (v(3:end,:)-2*v(2:end-1,:)+v(1:end-2,:))/dx^2;
vyy(:,2:end-1) = (v(:,3:end)-2*v(:,2:end-1)+v(:,1:end-2))/dy^2;

adv_u = up.*ux + vp.*uy;
adv_u = 0.5*(adv_u(:,1:Nj)+adv_u(:,2:Nj+1));

adv_v = up.*vx + vp.*vy;
adv_v = 0.5*(adv_v(1:Ni,:)+adv_v(2:Ni+1,:));

diff_u = nu*(uxx + uyy);
diff_v = nu*(vxx + vyy);

rhs_u = -0*adv_u + diff_u;
rhs_v = -0*adv_v + diff_v;

if(n==1)
    clear Cu Cv

    Cu{1} = rhs_u;
    Cv{1} = rhs_v;
    
    u = u + dt*Cu{1};
    v = v + dt*Cv{1};
elseif(n==2)
    Cu{2} = Cu{1};
    Cv{2} = Cv{1};
    Cu{1} = rhs_u;
    Cv{1} = rhs_v;

    u = u + dt*(3/2*Cu{1} - Cu{2});
    v = v + dt*(3/2*Cv{1} - Cv{2});
else
    Cu{3} = Cu{2};
    Cv{3} = Cv{2};
    Cu{2} = Cu{1};
    Cv{2} = Cv{1};
    Cu{1} = rhs_u;
    Cv{1} = rhs_v;

    u = u + dt*(23/12*Cu{1} - 4/3*Cu{2} + 5/12*Cu{3});
    v = v + dt*(23/12*Cv{1} - 4/3*Cv{2} + 5/12*Cv{3});
end    







