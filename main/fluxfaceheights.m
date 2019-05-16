function [Huwx,Huwy]=fluxfaceheights(H,d,h,u,v)

global inu_ij in_ij in_im1j inv_ij in_ijm1 out_im1j out_ij out_ijm1 ...
    outu_ij outv_ij outu_ip1j outv_ijp1

[Ni,Nj]=size(H);
Huwx=zeros(Ni+1,Nj);
Huwy=zeros(Ni,Nj+1);
  
Huwx(inu_ij)=0.5*(d(in_im1j)+d(in_ij))+(u(inu_ij)>0).*h(in_im1j)+(u(inu_ij)<0).*h(in_ij)+(u(inu_ij)==0).*max(h(in_im1j),h(in_ij));
Huwy(inv_ij)=0.5*(d(in_ijm1)+d(in_ij))+(v(inv_ij)>0).*h(in_ijm1)+(v(inv_ij)<0).*h(in_ij)+(v(inv_ij)==0).*max(h(in_ijm1),h(in_ij));

Huwx(Ni+1,:)=H(Ni,:);
Huwy(:,Nj+1)=H(:,Nj);

Huwx(outu_ij)=max(H(out_im1j),H(out_ij));
Huwy(outv_ij)=max(H(out_ijm1),H(out_ij));  

Huwx(outu_ip1j)=max(H(out_im1j),H(out_ij));
Huwy(outv_ijp1)=max(H(out_ijm1),H(out_ij));    
