function setup_pointers(Ni,Nj)

global_pointers;

getvariables;

% Pointers for fs2d vectorized
[i,j]=find(cellmark==0);
out_ij=sub2ind([Ni Nj],i,j);

out_ip1j=sub2ind(size(cellmark),min(i+1,Ni),j);
out_im1j=sub2ind(size(cellmark),max(i-1,1),j);
out_ijp1=sub2ind(size(cellmark),i,min(j+1,Nj));
out_ijm1=sub2ind(size(cellmark),i,max(j-1,1));

outu_ij=sub2ind([Ni+1 Nj],i,j);
outu_ip1j=sub2ind([Ni+1 Nj],i+1,j);
outu_im1j=sub2ind([Ni+1 Nj],max(i-1,1),j);

outv_ij=sub2ind([Ni Nj+1],i,j);
outv_ijp1=sub2ind([Ni Nj+1],i,j+1);
outv_ijm1=sub2ind([Ni Nj+1],i,max(j-1,1));

[i,j]=find(cellmark~=0);
in_ij=sub2ind(size(cellmark),i,j);

in_ip1j=sub2ind(size(cellmark),min(i+1,Ni),j);
in_im1j=sub2ind(size(cellmark),max(i-1,1),j);
in_ijp1=sub2ind(size(cellmark),i,min(j+1,Nj));
in_ijm1=sub2ind(size(cellmark),i,max(j-1,1));
% Not used
%in_ip1jm1=sub2ind(size(cellmark),min(i+1,Ni),max(j-1,1));

inu_ij=sub2ind([Ni+1 Nj],i,j);
inu_ip1j=sub2ind([Ni+1 Nj],min(i+1,Ni+1),j);
% Not used
%inu_im1j=sub2ind([Ni+1 Nj],max(i-1,1),j);
%inu_ijm1=sub2ind([Ni+1 Nj],i,max(j-1,1));
%inu_ip1jm1=sub2ind([Ni+1 Nj],min(i+1,Ni+1),max(j-1,1));

inv_ij=sub2ind([Ni Nj+1],i,j);
inv_ijp1=sub2ind([Ni Nj+1],i,min(j+1,Nj+1));
% Not used
%inv_ijm1=sub2ind([Ni Nj+1],i,max(j-1,1));
%inv_ip1jm1=sub2ind([Ni Nj+1],min(i+1,Ni),max(j-1,1));
%inv_ip1j=sub2ind([Ni Nj+1],min(i+1,Ni),j);

poly_ij = sub2ind([Ni+1 Nj+1],ipoly,jpoly);
edge_u_ij = sub2ind([Ni+1 Nj],ipoly,min(jpoly,Nj));
edge_u_ijm1 = sub2ind([Ni+1 Nj],ipoly,max(jpoly-1,1));
edge_v_ij = sub2ind([Ni Nj+1],min(ipoly,Ni),jpoly);
edge_v_im1j = sub2ind([Ni Nj+1],max(ipoly-1,1),jpoly);

if((exist('ADVECTION','var') & ADVECTION) | ...
   (exist('DIFFUSION') & DIFFUSION))
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Pointers for u-face advection/diffusion
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [is,js]=find(markuface==0);
    advu_ij = sub2ind([Ni+1 Nj],is,js);
    bcu_ip1j = ones(size(advu_ij));
    bcu_im1j = ones(size(advu_ij));
    bcu_ijp1 = ones(size(advu_ij));
    bcu_ijm1 = ones(size(advu_ij));
    im1j = is-1;
    ip1j = is+1;
    ijm1 = js-1;
    ijp1 = js+1;
    
    % bc_grad: 0 = free slip, 1 = no slip
    if(~exist('NOSLIP','var') | ~NOSLIP)
        bc_grad = 0;
    else
        bc_grad = 1; 
    end
    
    % These are already taken care of b/c advu_ip1j<=Ni+1 and
    % advu_im1j>=1 since 2<=advu_ij<=Ni, i.e. markuface does not
    % include boundary edges.
    % dudx(advu_ij) = (bcu_ip1j.*u(advu_ip1j)-bcu_im1j.*u(advu_im1j))/(2*dx);
    
    % These need to be taken care of 
    % dudy(advu_ij) = (bcu_ijp1.*u(advu_ijp1)-bcu_ijm1.*u(advu_ijm1))/(2*dy);
    for m=1:length(is)
        if(js(m)==1)
            ijm1(m) = 1;
            if(bc_grad==0) % Free slip
                bcu_ijm1(m) = 1;
            else % No slip
                bcu_ijm1(m) = -1;
            end
        elseif(js(m)==Nj)
            ijp1(m) = Nj;
            if(bc_grad==0) % Free slip
                bcu_ijp1(m) = 1;
            else % No slipe
                bcu_ijp1(m) = -1;
            end
        elseif(cellmark(is(m),js(m)-1)==0)
            ijm1(m) = js(m);
            if(bc_grad==0) % Free slip
                bcu_ijm1(m) = 1;
            else
                bcu_ijm1(m) = -1;
            end
        elseif(cellmark(is(m),js(m)+1)==0)
            ijp1(m) = js(m);
            if(bc_grad==0) % Free slip
                bcu_ijp1(m) = 1;
            else
                bcu_ijp1(m) = -1;
            end
        end
    end
    advu_ip1j = sub2ind([Ni+1 Nj],ip1j,js);
    advu_im1j = sub2ind([Ni+1 Nj],im1j,js);
    advu_ijp1 = sub2ind([Ni+1 Nj],is,ijp1);
    advu_ijm1 = sub2ind([Ni+1 Nj],is,ijm1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Pointers for v-face advection/diffusion
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [is,js]=find(markvface==0);
    advv_ij = sub2ind([Ni Nj+1],is,js);
    bcv_ip1j = ones(size(advv_ij));
    bcv_im1j = ones(size(advv_ij));
    bcv_ijp1 = ones(size(advv_ij));
    bcv_ijm1 = ones(size(advv_ij));
    im1j = is-1;
    ip1j = is+1;
    ijm1 = js-1;
    ijp1 = js+1;
    
    % These are already taken care of b/c advv_ijm1>=1 and
    % advv_ijp1<=Nj since 2<=advv_ij<=Nj, i.e. markvface does not
    % include boundary edges.
    % dvdy(advv_ij) = (bcv_ijp1.*v(advv_ijp1)-bcv_ijm1.*v(advv_ijm1))/(2*dy);
    
    % These need to be taken care of 
    % dvdx(advv_ij) = (bcv_ip1j.*v(advv_ip1j)-bcv_im1j.*v(advv_im1j))/(2*dx);
    for m=1:length(is)
        if(is(m)==1)
            im1j(m) = 1;
            if(bc_grad==0) % Free slip
                bcv_im1j(m) = 1;
            else % No slip
                bcv_im1j(m) = -1;
            end
        elseif(is(m)==Ni)
            ip1j(m) = Ni;
            if(bc_grad==0) % Free slip
                bcv_ip1j(m) = 1;
            else % No slipe
                bcv_ip1j(m) = -1;
            end
        elseif(cellmark(is(m)-1,js(m))==0)
            im1j(m) = is(m);
            if(bc_grad==0) % Free slip
                bcv_im1j(m) = 1;
            else
                bcv_im1j(m) = -1;
            end
        elseif(cellmark(is(m)+1,js(m))==0)
            ip1j(m) = is(m);
            if(bc_grad==0) % Free slip
                bcv_ip1j(m) = 1;
            else
                bcv_ip1j(m) = -1;
            end
        end
    end
    advv_ip1j = sub2ind([Ni Nj+1],ip1j,js);
    advv_im1j = sub2ind([Ni Nj+1],im1j,js);
    advv_ijp1 = sub2ind([Ni Nj+1],is,ijp1);
    advv_ijm1 = sub2ind([Ni Nj+1],is,ijm1);

    % Need to set the u and v velocities on the type 3 edge to that of
    % the interior edge, since the type 3 edge is not computed
    %u(u_type3)=u(u_type3_inner);
    %v(v_type3)=v(v_type3_inner);
    
    [is,js]=find(markuface==3);
    iu1=zeros(size(is));
    iu2=zeros(size(is));
    for m=1:length(is)
        if(is(m)==1)
            iu1(m)=2;
            iu2(m)=3;
        elseif(is(m)==Ni+1)
            iu1(m)=Ni;
            iu2(m)=Ni-1;
        elseif(cellmark(is(m)-1,js(m))==0)
            iu1(m)=is(m)+1;
            iu2(m)=is(m)+2;
        else
            iu1(m)=is(m);
            iu2(m)=is(m)-1;
        end
    end
    u_type3u = sub2ind([Ni+1 Nj],is,js);
    v_type3u = sub2ind([Ni Nj+1],min(is,Ni),js);
    
    u_type3u_inner1=sub2ind([Ni+1 Nj],iu1,js);
    u_type3u_inner2=sub2ind([Ni+1 Nj],iu2,js);
    
    v_type3u_inner1=sub2ind([Ni Nj+1],min(iu1,Ni),js);
    v_type3u_inner2=sub2ind([Ni Nj+1],min(iu2,Ni),js);
    
    v_type3u(markvface(v_type3u_inner1)~=0)=-1;
    v_type3u(markvface(v_type3u_inner2)~=0)=-1;
    
    v_type3u_inner1(v_type3u==-1)=[];
    v_type3u_inner2(v_type3u==-1)=[];
    v_type3u(v_type3u==-1)=[];
    
    [is,js]=find(markvface==3);
    jv1=zeros(size(is));
    jv2=zeros(size(is));
    for m=1:length(is)
        if(js(m)==1)
            jv1(m)=2;
            jv2(m)=3;
        elseif(js(m)==Nj+1)
            jv1(m)=Nj;        
            jv2(m)=Nj-1;
        elseif(cellmark(is(m),js(m)-1)==0)
            jv1(m)=js(m)+1;
            jv2(m)=js(m)+2;
        else
            jv1(m)=js(m);
            jv2(m)=js(m)-1;
        end
    end
    u_type3v = sub2ind([Ni+1 Nj],is,min(js,Nj));
    v_type3v = sub2ind([Ni Nj+1],is,js);
    
    v_type3v_inner1=sub2ind([Ni Nj+1],is,jv1);
    v_type3v_inner2=sub2ind([Ni Nj+1],is,jv2);
    
    u_type3v_inner1=sub2ind([Ni+1 Nj],is,min(jv1,Nj));
    u_type3v_inner2=sub2ind([Ni+1 Nj],is,min(jv2,Nj));
    
    u_type3v(markuface(u_type3v_inner1)~=0)=-1;
    u_type3v(markuface(u_type3v_inner2)~=0)=-1;
    
    u_type3v_inner1(u_type3v==-1)=[];
    u_type3v_inner2(u_type3v==-1)=[];
    u_type3v(u_type3v==-1)=[];
end


