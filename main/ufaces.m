function [u_vface,v_uface]=ufaces(u,v)

global markuface markvface

[Nip1,Nj]=size(u);
Ni=Nip1-1;
u_vface=zeros(Ni,Nj+1);
v_uface=zeros(Ni+1,Nj);

u_vface(:,2:Nj-1)=0.25*(u(1:Ni,1:Nj-2)+u(1:Ni,2:Nj-1)+u(2:Ni+1,1:Nj-2)+u(2:Ni+1,2:Nj-1));
v_uface(2:Ni-1,:)=0.25*(v(1:Ni-2,1:Nj)+v(2:Ni-1,1:Nj)+v(1:Ni-2,2:Nj+1)+v(2:Ni-1,2:Nj+1));

if 0
[is,js]=find(markuface~=0);
for m=1:length(is)
    i00=is(m)+[1,0,-1,0];
    j00=is(m)+[0,1,0,-1];
    i0=i00;
    j0=j00;
    i00(i0<1 | i0>Ni | j0<1 | j0>Nj+1)=[];
    j00(i0<1 | i0>Ni | j0<1 | j0>Nj+1)=[];
    i0=i00;
    j0=j00;
    v_uface(is(m),js(m))=max(max(v(i0,j0)));
end

[is,js]=find(markvface~=0);
for m=1:length(is)
    i00=is(m)+[1,0,-1,0];
    j00=is(m)+[0,1,0,-1];
    i0=i00;
    j0=j00;
    i00(i0<1 | i0>Ni+1 | j0<1 | j0>Nj)=[];
    j00(i0<1 | i0>Ni+1 | j0<1 | j0>Nj)=[];
    i0=i00;
    j0=j00;
    u_vface(is(m),js(m))=max(max(u(i0,j0)));
end
end

%v_uface(inu_ij)=0.25*(v(inv_ijm1)+v(inv_ij)+v(inv_ip1jm1)+v(inv_ip1j));
%u_vface(inv_ij)=0.25*(u(inu_ijm1)+u(inu_ij)+u(inu_ip1jm1)+u(inu_ip1j));
