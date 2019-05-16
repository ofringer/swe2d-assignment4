% Cell-centered velocities
uc = 0.5*(u(1:end-1,:)+u(2:end,:));
vc = 0.5*(v(:,1:end-1)+v(:,2:end));

if(n==1)
    nouts=[nout:nout:nsteps];
    if(~isempty(find('h'==vars_to_output)))
        h_out=zeros(length(nouts),length(x_transect));
    end
    if(~isempty(find('u'==vars_to_output)))    
        u_out=zeros(length(nouts),length(x_transect));
    end
    if(~isempty(find('v'==vars_to_output)))    
        v_out=zeros(length(nouts),length(x_transect));
    end
    if(~isempty(find('C'==vars_to_output)))        
        C_out=zeros(length(nouts),length(x_transect));
    end
    if(~isempty(find('H'==vars_to_output)))        
        H_out=zeros(length(nouts),length(x_transect));
    end
    if(~isempty(find('z'==vars_to_output)))        
        zb_out=zeros(length(nouts),length(x_transect));
    end
end

if(nout==1 | (rem(n,nout)==0 & n>1))
    kout=find(n==nouts);
    for m=1:length(vars_to_output)
        var_p=vars_to_output(m);
        switch(var_p)
          case 'h'
            h_out(kout,:)=interp2(xc',yc',h',x_transect,y_transect)';
          case 'u'
            u_out(kout,:)=interp2(xc',yc',uc',x_transect,y_transect)';
          case 'v'
            v_out(kout,:)=interp2(xc',yc',vc',x_transect,y_transect)';
          case 'C'
            C_out(kout,:)=interp2(xc',yc',C',x_transect,y_transect)';
          case 'H'
            H_out(kout,:)=interp2(xc',yc',h'+d',x_transect,y_transect)';
          case 'z'
            zb_out(kout,:)=interp2(xc',yc',zb',x_transect,y_transect)';
          otherwise
            error(sprintf('No variable type %s to output.',var_p));
        end
    end
end