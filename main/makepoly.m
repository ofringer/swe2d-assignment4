function [xpoly,ypoly,ipoly,jpoly,inp_ij]=makepoly(dx,dy,cellmark)

[Ni,Nj]=size(cellmark);
[xp,yp]=ndgrid([0:Ni]*dx,[0:Nj]*dy);

p = false(Ni+1,Nj+1);
inp_ij = false(Ni+1,Nj+1);
for i=1:Ni
    for j=1:Nj
        if(cellmark(i,j)~=0)
            inp_ij(i,j)=true;
            inp_ij(i+1,j)=true;
            inp_ij(i,j+1)=true;

            if(i==1 | cellmark(i-1,j)==0)
                p(i,j)=true;
                p(i,j+1)=true;
            end
            if(i==Ni | cellmark(i+1,j)==0)
                p(i+1,j)=true;
                p(i+1,j+1)=true;
            end
            if(j==1 | cellmark(i,j-1)==0)
                p(i,j)=true;
                p(i+1,j)=true;
            end
            if(j==Nj | cellmark(i,j+1)==0)
                p(i,j+1)=true;
                p(i+1,j+1)=true;
            end
        end
    end
end

xp0=xp;
yp0=yp;
xp=xp(find(p==true));
yp=yp(find(p==true));
[ip,jp]=find(p==true);

ip=ip(:);
jp=jp(:);
xp=xp(:);
yp=yp(:);
N=length(ip);
found=false(N,1);

all_ind=1;
found(all_ind)=true;
idir=[1,0,-1,0];
jdir=[0,1,0,-1];
for n=1:N
    for m=1:4
        ind=find(ip==ip(all_ind(end))+idir(m) & ...
                 jp==jp(all_ind(end))+jdir(m));
        if(~found(ind))
            all_ind = [all_ind,ind];
            found(ind) = true;
            break;
        end
    end
end

xpoly=xp(all_ind);
ypoly=yp(all_ind);
ipoly=ip(all_ind);
jpoly=jp(all_ind);





