function [uface,vface]=faces_to_zero(uface,vface,outu_ij,outv_ij,outu_ip1j,outv_ijp1,markuface,markvface)

uf0=uface;
vf0=vface;

[Nip1,Nj]=size(uface);
Ni=Nip1-1;

uface(Ni+1,:)=0;
vface(:,Nj+1)=0;

uface(outu_ij)=0;
vface(outv_ij)=0;

uface(outu_ip1j)=0;
vface(outv_ijp1)=0;

uface(markuface~=0)=0;
vface(markvface~=0)=0;
uface(markuface==2)=uf0(markuface==2);
vface(markvface==2)=vf0(markvface==2);
