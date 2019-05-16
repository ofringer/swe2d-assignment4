function phi = boundaryconditions_cell(phi,xc,yc,ntime,boundaryregion,phi_b)

global cellmark

ntime0=ntime;
for m=1:length(phi_b)
  bound=boundaryregion{m};
  if(~isempty(bound))
    phi_indices=find(cellmark==3 & inpolygon(xc,yc,bound(:,1),bound(:,2)));
  else
    phi_indices=find(cellmark==3);
  end

  if(length(phi_b{m})==1)
      ntime=1;
  else
      ntime=ntime0;
  end
  phi(phi_indices)=phi_b{m}(ntime);
end


  





