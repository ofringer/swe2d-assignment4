function progress(u,v,C,t_start,n,dt,dx,dy,nsteps,Cmax_allowed,C0_max_allowed)

  global Cmax0

  % Maximum Courant number
  uc=0.5*(u(1:end-1,:)+u(2:end,:));
  vc=0.5*(v(:,1:end-1)+v(:,2:end));
  Cmax=max(max((abs(uc)*dt/dx+abs(vc)*dt/dy)));

  if(n==1)
    Cmax0=0;
  else
    if(Cmax>Cmax0)
      Cmax0=Cmax;
    end
  end

  if(Cmax>Cmax_allowed)
      error(sprintf('Time step: %d, Maximum allowable Courant number = %.2f>%.2f.',n,Cmax,Cmax_allowed));
  end
  if(~isempty(C0_max_allowed) & max(max(C))>C0_max_allowed)
      error(sprintf('Time step: %d, Maximum allowable SSC = %.2f>%.2f.',n,max(max(C)),C0_max_allowed));
  end
  if(~isempty(find(isnan(C))) | ~isempty(find(isnan(u))) | ~isempty(find(isnan(v))))
      error(sprintf('Time step: %d, Model is blowing up!\n',n));
  end

  elapsed_time=toc(t_start);
  time_per_step=elapsed_time/n;
  time_remaining=(nsteps-n)*time_per_step;

  fprintf('\n');
  fprintf('Time step: %d of %d (Cmax=%.2f), Simulation time: %.2f hr\n',n,nsteps,Cmax,n*dt/3600);
  fprintf('Wallclock time per time step: %.2f s, Speedup: %.2f X real time\n', ...
	  time_per_step,dt/time_per_step);
  fprintf('Time remaining: %.0f s (%.1f min)\n',time_remaining,time_remaining/60);
  if(n==nsteps)
      fprintf('\n\nTotal Wallclock time: %.2f s, %.2e s per step.\n',...
          elapsed_time,elapsed_time/nsteps);
      end
  fprintf('\n\n\n\n');
