global C_vals Q_vals A_vals U_vals

getvariables;

load(flux_file);

if 0, % For debugging
xbox=zeros(Ni+1,Nj+1);
ybox=zeros(Ni+1,Nj+1);
cbox=zeros(Ni+1,Nj+1);

xbox(1:Ni,1:Nj)=xc-dx/2;
ybox(1:Ni,1:Nj)=yc-dy/2;
xbox(end,1:Nj)=xbox(end-1,1:Nj)+dx;
ybox(end,1:Nj)=ybox(end-1,1:Nj);
xbox(1:Ni+1,Nj+1)=xbox(1:Ni+1,Nj);
ybox(1:Ni+1,Nj+1)=ybox(1:Ni+1,Nj)+dy;
%[xbox,ybox]=ndgrid([1:Ni+1],[1:Nj+1]);
uc(cellmark==0)=nan;
Cbox(1:Ni,1:Nj)=uc;
Cbox(Ni+1,1:Nj)=uc(Ni,1:Nj);
Cbox(:,Nj+1)=Cbox(:,Nj);

figure(1)
clf;
hold on;
pcolor(xbox,ybox,Cbox);

for j=1:5
    if(flux_type(j)=='u')
        xs=xu(flux_i{j},flux_j{j});
        ys=yu(flux_i{j},flux_j{j});
    else
        xs=xv(flux_i{j},flux_j{j});
        ys=yv(flux_i{j},flux_j{j});
    end
    plot(xs,ys,'r.','markersize',10);
end
end

% Compute fluxes here
if(n==1)
    for mf=1:length(flux_i)
        C_vals{mf}=zeros(nsteps,max([length(flux_i{mf}),length(flux_j{mf})]));
        Q_vals{mf}=zeros(nsteps,max([length(flux_i{mf}),length(flux_j{mf})]));
        U_vals{mf}=zeros(nsteps,max([length(flux_i{mf}),length(flux_j{mf})]));
        A_vals{mf}=zeros(nsteps,max([length(flux_i{mf}),length(flux_j{mf})]));
    end
end

[Ni,Nj]=size(C);
for mf=1:length(flux_i)
    is=flux_i{mf};
    js=flux_j{mf};
    inds_u = sub2ind([Ni+1 Nj],is,min(js,Nj));
    inds_v = sub2ind([Ni Nj+1],min(is,Ni),js);
    C_vals{mf}(n,:)=(flux_type{mf}=='u').*Cuwx(inds_u)+(flux_type{mf}=='v').*Cuwy(inds_v);
    Q_vals{mf}(n,:)=-(flux_type{mf}=='u').*Huwx(inds_u).*utheta(inds_u)*dy ...
        + (flux_type{mf}=='v').*Huwy(inds_v).*vtheta(inds_v)*dx;
    U_vals{mf}(n,:)=-(flux_type{mf}=='u').*utheta(inds_u)+(flux_type{mf}=='v').*vtheta(inds_v);
    A_vals{mf}(n,:)=(flux_type{mf}=='u').*Huwx(inds_u).*dy+(flux_type{mf}=='v').*Huwy(inds_v).*dx;
end
