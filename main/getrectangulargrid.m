function [dx,dy,Ni,Nj,depth,markuface,markvface,cellmark]=getrectangulargrid()

Ni=100;
Nj=20;
L=100e3;
W=20e3;
dx=L/Ni;
dy=W/Nj;
depth=10*ones(Ni,Nj);

Lstep=30e3;
Hstep=10e3;

[xc,yc]=ndgrid([dx/2:dx:L-dx/2],[dy/2:dy:W-dy/2]);
% 0 = outside
% 1 = inside
% 3 = type 3 free-surface forced
cellmark=ones(Ni,Nj);
cellmark(end,:)=3;
cellmark(xc<=Lstep & yc<=Hstep)=0;

% 0 computational edge
% 1 no-flux edge
% 2 flow-specified edge
markuface=zeros(Ni+1,Nj);
markvface=zeros(Ni,Nj+1);

markuface(1,:)=2;
markuface(end,:)=3;

[i,j]=find(cellmark==0);
out_ij=sub2ind(size(cellmark),i,j);

outu_ij=sub2ind([Ni+1 Nj],i,j);
outu_ip1j=sub2ind([Ni+1 Nj],i+1,j);
outu_im1j=sub2ind([Ni+1 Nj],max(i-1,1),j);

outv_ij=sub2ind([Ni Nj+1],i,j);
outv_ijp1=sub2ind([Ni Nj+1],i,j+1);
outv_ijm1=sub2ind([Ni Nj+1],i,max(j-1,1));

markuface(outu_ij)=0;
markuface(outu_ip1j)=0;
markvface(outv_ij)=0;
markvface(outv_ijp1)=0;

