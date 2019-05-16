function plotswe(xc,yc,u,v,h,C,zb,d,cellmark,n,nsteps,dt,vars_to_plot,MOVIE)

fontsize=12;

global swe_movie n_movie

[Ni,Nj]=size(C);
% Show log of scalar-mean(scalar) when (max(scalar)-min(scalar))/mean(scalar)<tol_for_plot
tol_for_plot=1e-4;

% Need L and W for axes
dx=xc(2,1)-xc(1,1);
dy=yc(1,2)-yc(1,1);
L=max(max(xc))+dx/2;
W=max(max(yc))+dy/2;

% Cell-centered velocities
uc = 0.5*(u(1:end-1,:)+u(2:end,:));
vc = 0.5*(v(:,1:end-1)+v(:,2:end));

hp=h;
up=uc;
vp=vc;
Cp=C;
zp=zb;
Hp=h+d;

hp(cellmark==0)=nan;
up(cellmark==0)=nan;
vp(cellmark==0)=nan;
Cp(cellmark==0)=nan;
zp(cellmark==0)=nan;
Hp(cellmark==0)=nan;

for m=1:length(vars_to_plot)
    var_p=vars_to_plot(m);
    switch(var_p)
      case 'h'
          var_plot{m}=hp;
          plot_title{m}='Free surface (m)';
      case 'H'
          var_plot{m}=Hp;
          plot_title{m}='Water depth (h+d; m)';
      case 'u'
        var_plot{m}=up;
        plot_title{m}='u velocity (m/s)';
      case 'v'
        var_plot{m}=vp;
        plot_title{m}='v velocity (m/s)';
      case 'C'
        var_plot{m}=Cp;
        plot_title{m}='SSC (kg/m^3)';
      case 'z'
        var_plot{m}=zp;
        plot_title{m}='Bed height (m)';
      case 'q'
        var_plot{m}=[];
        plot_title{m}='Vectors';
      case 's'
        var_plot{m}=[];
        plot_title{m}='Streamlines';
      otherwise
        error(sprintf('No variable type %s to plot.',var_p));
    end
end

figure(1);
clf;
hold on;
set(gca,'fontsize',fontsize,'box','on');
for m=1:length(vars_to_plot)
    var_p=vars_to_plot(m);
    if(length(vars_to_plot)>3)
        subplot(3,2,m);
    else
        subplot(length(vars_to_plot),1,m);
    end
    hold on;
    set(gca,'box','on')
    if(var_p=='q')
        quiver(xc,yc,up,vp,2,'k-');
    elseif(var_p=='s')
        iskip=fix(Ni/20);
        jskip=fix(Nj/20);
        xs=xc([iskip:iskip:Ni],[jskip:jskip:Nj]);
        ys=yc([iskip:iskip:Ni],[jskip:jskip:Nj]);
        hs=streamline(stream2(xc',yc',uc',vc',xs(:),ys(:),[0.1 100]));
        for ms=1:length(hs)
            hs(ms).Color='k';
        end
    else
        vplot=var_plot{m};
        vplot=vplot(~isnan(vplot));
        vplot=vplot(:);
        if(abs((max(vplot)-min(vplot))/max(vplot))<tol_for_plot)
            phi_plot=log(abs(var_plot{m}-mean(vplot)));
        else
            phi_plot=var_plot{m};
        end
        pcolor(xc,yc,phi_plot);
    end
    title(sprintf('%s, Step %d of %d, time=%.2f hr',plot_title{m},n,nsteps,n*dt/3600),'fontsize',fontsize);
    shading flat;
    if(var_p~='q' & var_p~='s')
        colorbar;
    end

    axis image;
    axis([0 L 0 W]);
    xlabel('x (m)')
    ylabel('y (m)')
    %title(sprintf('Free surface (m) and vectors at time-step %d
    %(t=%.1f s)\n',n,t(n)));
end
drawnow;

if(MOVIE)
    figure(gcf);
    if(n==1)
        n_movie=1;
        swe_movie=getframe;
    else
        swe_movie(n_movie)=getframe(gcf);
    end
    n_movie=n_movie+1;
end

if(0) % Line plot for debugging
    [Ni,Nj]=size(h);
    x=xc(:,fix(Nj/2));
    figure(2)
    subplot(3,1,1)
    Hi = h(:,fix(Nj/2))+d(:,fix(Nj/2));
    plot(x,h(:,4));
    subplot(3,1,2)
    plot([1:Ni+1],u(:,[1:10]))
    subplot(3,1,3)
    plot([1:Ni],v(:,[1:10]));
    drawnow;
end