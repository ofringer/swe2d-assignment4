%
% Set the velocity at the boundaries given the inflow fluxes.
%
% Oliver Fringer
% Yun Zhang
% Stanford University
%
function [flux_x,flux_y]=...
boundaryconditions_face(phi,u,v,flux_x,flux_y,A_uface,A_vface,xu,yu,xv,yv,ntime,boundaryregion,Qb,ub,vb,phi_b,type)

global markuface markvface u_indices v_indices

[Ni,Nj]=size(flux_x);
Ni=Ni-1;

for m=1:length(boundaryregion)
  bound=boundaryregion{m};
  if(~isempty(bound))
    u_indices=find(markuface==2 & inpolygon(xu,yu,bound(:,1),bound(:,2)));
    v_indices=find(markvface==2 & inpolygon(xv,yv,bound(:,1),bound(:,2)));
  else
    u_indices=find(markuface==2);
    v_indices=find(markvface==2);
  end

  if(type=='Q')
      A=sum(A_uface(u_indices))+sum(A_vface(v_indices));
      if(length(Qb{m})==1)
          ntime=1;
      end
      Uave=Qb{m}(ntime)/A;
      
      flux_x(u_indices)=Uave;
      flux_y(v_indices)=Uave;
  elseif(type=='u')
      if(length(ub{m})==1)
          ntime=1;
      end

      flux_x(u_indices)=ub{m}(ntime);
      flux_y(v_indices)=vb{m}(ntime);
  else
      if(length(phi_b{m})==1)
          ntime=1;
      end
      [ix,jx]=ndgrid(1:Ni+1,1:Nj);
      [iy,jy]=ndgrid(1:Ni,1:Nj+1);

      ix_plus = ix(u_indices(u(u_indices)>=0));
      jx_plus = jx(u_indices(u(u_indices)>=0));
      ic=ix_plus-1;
      jc=jx_plus;
      for k=1:length(ic)
          if(ic(k)<=0 | ic(k)>Ni)
              flux_x(ix_plus(k),jx_plus(k))=phi_b{m}(ntime);
          else
              flux_x(ix_plus(k),jx_plus(k))=phi(ic(k),jc(k));
          end
      end

      ix_minus = ix(u_indices(u(u_indices)<0));
      jx_minus = jx(u_indices(u(u_indices)<0));
      ic=ix_minus;
      jc=jx_minus;
      for k=1:length(ic)
          if(ic(k)<=0 | ic(k)>Ni)
              flux_x(ix_minus(k),jx_minus(k))=phi_b{m}(ntime);
          else
              flux_x(ix_minus(k),jx_minus(k))=phi(ic(k),jc(k));
          end
      end

      iy_plus = iy(v_indices(v(v_indices)>=0));
      jy_plus = jy(v_indices(v(v_indices)>=0));
      ic=iy_plus;
      jc=jy_plus-1;
      for k=1:length(ic)
          if(jc(k)<=0 | jc(k)>Nj)
              flux_y(iy_plus(k),jy_plus(k))=phi_b{m}(ntime);
          else
              flux_y(iy_plus(k),jy_plus(k))=phi(ic(k),jc(k));
          end
      end

      iy_minus = iy(v_indices(v(v_indices)<0));
      jy_minus = jy(v_indices(v(v_indices)<0));
      ic=iy_minus;
      jc=jy_minus;
      for k=1:length(ic)
          if(jc(k)<=0 | jc(k)>Nj)
              flux_y(iy_minus(k),jy_minus(k))=phi_b{m}(ntime);
          else
              flux_y(iy_minus(k),jy_minus(k))=phi(ic(k),jc(k));
          end
      end

      %      flux_x(u_indices)=phi_b{m}(ntime);
      %      flux_y(v_indices)=phi_b{m}(ntime);
  end      

  if(type=='Q')
      % Set boundary velocities to negative if flow inward from top or right
      for m=1:length(u_indices)
          [i,j]=ind2sub(size(flux_x),u_indices(m));
          if(i==length(flux_x(:,1)) | (i<=Ni & cellmark(i,j)==0))
              flux_x(i,j)=-flux_x(i,j);
          end
      end
      
      for m=1:length(v_indices)
          [i,j]=ind2sub(size(flux_y),v_indices(m));
          if(j==length(flux_y(1,:)) | (j<=Nj & cellmark(i,j)==0))
              flux_y(i,j)=-flux_y(i,j);
          end
      end
  end
end
  





