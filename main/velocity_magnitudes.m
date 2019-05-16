%
% Compute the magnitudes of the U and V velocities on the edges for
% the drag term.
%
% Oliver Fringer
% Yun Zhang
% Stanford University
%
function [Ux,Uy]=velocity_magnitudes(u,v,h,Huwx,Huwy,g,dx,dy,dt,pointers)

getpointers;

[Ni,Nj]=size(h);
utilde=zeros(Ni+1,Nj);
vtilde=zeros(Ni,Nj+1);

% Provisional value of u for use in the drag term inside the domain
utilde(inu_ij)=u(inu_ij)-g*dt/(2*dx)*(h(in_ij)-h(in_im1j));
vtilde(inv_ij)=v(inu_ij)-g*dt/(2*dy)*(h(in_ij)-h(in_ijm1));

[utilde,vtilde]=faces_to_zero(utilde,vtilde,outu_ij,outv_ij,outu_ip1j,outv_ijp1,markuface,markvface);
[u_vface,v_uface]=ufaces(utilde,vtilde);
faces_to_zero(u_vface,v_uface,outu_ij,outv_ij,outu_ip1j,outv_ijp1,markuface,markvface);

Ux=sqrt(u.^2+v_uface.^2);
Uy=sqrt(u_vface.^2+v.^2);
