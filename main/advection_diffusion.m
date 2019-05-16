function [u,v,Cu,Cv] = advection_diffusion(u,v,dx,dy,Cu,Cv,n,dt,nu,ADVECTION,DIFFUSION)

global_pointers;

[Ni,Nj]=size(cellmark);

rhs_u = zeros(Ni+1,Nj);
rhs_v = zeros(Ni,Nj+1);

if(ADVECTION)
    adv_u = zeros(Ni+1,Nj);
    adv_v = zeros(Ni,Nj+1);
    
    v_uface = zeros(Ni+1,Nj);
    v_uface(2:Ni,:) = 0.25*(v(1:Ni-1,1:Nj)+v(2:Ni,1:Nj)+v(1:Ni-1,2:Nj+1)+v(2:Ni,2:Nj+1));
    
    u_vface = zeros(Ni,Ni+1);
    u_vface(:,2:Nj) = 0.25*(u(1:Ni,1:Nj-1)+u(2:Ni+1,1:Nj-1)+u(1:Ni,2:Nj)+u(2:Ni+1,2:Nj));
    
    adv_u(advu_ij) = u(advu_ij).*(bcu_ip1j.*u(advu_ip1j)-bcu_im1j.*u(advu_im1j))/(2*dx) + ...
        v_uface(advu_ij).*(bcu_ijp1.*u(advu_ijp1)-bcu_ijm1.*u(advu_ijm1))/(2*dy);

    adv_v(advv_ij) = u_vface(advv_ij).*(bcv_ip1j.*v(advv_ip1j)-bcv_im1j.*v(advv_im1j))/(2*dx) + ...
        v(advv_ij).*(bcv_ijp1.*v(advv_ijp1)-bcv_ijm1.*v(advv_ijm1))/(2*dy);

    rhs_u = rhs_u - adv_u;
    rhs_v = rhs_v - adv_v;
end

if(DIFFUSION)
    diff_u = zeros(Ni+1,Nj);
    diff_v = zeros(Ni,Nj+1);

    diff_u(advu_ij) = nu*(bcu_ip1j.*u(advu_ip1j)-2*u(advu_ij)+bcu_im1j.*u(advu_im1j))/dx^2 + ...
        nu*(bcu_ijp1.*u(advu_ijp1)-2*u(advu_ij)+bcu_ijm1.*u(advu_ijm1))/dy^2;
    diff_v(advv_ij) = nu*(bcv_ip1j.*v(advv_ip1j)-2*v(advv_ij)+bcv_im1j.*v(advv_im1j))/dx^2 + ...
        nu*(bcv_ijp1.*v(advv_ijp1)-2*v(advv_ij)+bcv_ijm1.*v(advv_ijm1))/dy^2;

    rhs_u = rhs_u + diff_u;
    rhs_v = rhs_v + diff_v;
end

rhs_u(u_type3u)=0;
rhs_u(u_type3u_inner1)=0;
rhs_v(v_type3u)=0;
%rhs_v(v_type3u_inner1)=0;
rhs_u(u_type3v)=0;
rhs_u(u_type3v_inner1)=0;
%rhs_v(v_type3v)=0;
rhs_v(v_type3v_inner1)=0;

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







