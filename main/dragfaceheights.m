function [Hdragx,Hdragy]=dragfaceheights(H)

global inu_ij in_im1j in_ij inv_ij in_ijm1 outu_ij out_im1j out_ij ...
    outv_ij out_ijm1 outu_ip1j outv_ijp1

[Ni,Nj]=size(H);
Hdragx=zeros(Ni+1,Nj);
Hdragy=zeros(Ni,Nj+1);
  
Hdragx(inu_ij)=0.5*(H(in_im1j)+H(in_ij));
Hdragy(inv_ij)=0.5*(H(in_ijm1)+H(in_ij));

Hdragx(Ni+1,:)=H(Ni,:);
Hdragy(:,Nj+1)=H(:,Nj);

Hdragx(outu_ij)=max(H(out_im1j),H(out_ij));
Hdragy(outv_ij)=max(H(out_ijm1),H(out_ij));  

Hdragx(outu_ip1j)=max(H(out_im1j),H(out_ij));
Hdragy(outv_ijp1)=max(H(out_ijm1),H(out_ij));    
