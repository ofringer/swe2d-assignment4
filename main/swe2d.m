%
% SWE2D 2D Shallow water equation solver
% 
% Grid size: Ni by Nj
%
% Edge markers: 0=computational, 1=closed, -1=inactive, 2=flux specified
%   markuface: Ni+1 by Nj markers for u-face edges
%   markvface: Ni by Nj+1 markers for v-face edges
%
% Cell marker: 0=inactive, 1=computational, 2=not used,
% 3=water-level specified
%   cellmark: Ni by Nj markers for cell-centered quanties -
%   scalars, free surface.
%
% Main input file: getvariables.m
% Directory structure: 
%   ROOT/main - contains all shallow water code
%   ROOT/examples/example_name - contains codes for initial and
%     boundary conditions. The code should be run with 
%     ROOT/examples/example_name/run_swe.m
%
% Module variables (set to true/false - default in parantheses):
%  SSC (true): Suspended sediment transport
%  EXNER (false): Solve exner equation
%  BEDLOAD (false): Bed-load transport
%  EROSION (false): Include erosion from the bed
%  DEPOSITION (false): Include deposition
%  ADVECTION (false): Eulerian momentum advection
%  DIFFUSION (false): Momentum diffusion
%  HELICAL (false): Helical flow model
%  HELICAL_ADVECTION (false): Advetion of helicity
%  DRAGPATCH (false): Variable bottom drag (for vegetation model)
%  

%
% Oliver Fringer
% Yun Zhang
% Stanford University
% 

%
% Contains all global variables needed for the pointers
%
global_pointers;

%
% Variables needed for model run, like time step size, # grid
% points, etc...
%
getvariables;

% 
% Set the boundary conditions here
%
boundary_conditions;

%
% Grid file must contain:
% dx, dy, Ni, Nj, depth, markuface, markvface, cellmark
%
[dx,dy,Ni,Nj,depth,markuface,markvface,cellmark]=make_grid();
[xpoly,ypoly,ipoly,jpoly,inp_ij]=makepoly(dx,dy,cellmark);

% 
% Set up the grid - Note that Ni, Nj, dx, and dy are obtained when
% creating the grid above.
%
L=Ni*dx;
W=Nj*dy;
[xu,yu]=ndgrid([0:dx:L],[dy/2:dy:W-dy/2]);
[xv,yv]=ndgrid([dx/2:dx:L-dx/2],[0:dy:W]);
[xc,yc]=ndgrid([dx/2:dx:L-dx/2],[dy/2:dy:W-dy/2]);

%
% Set the markuface and markvface markers based on the boundary
% type
%
boundary_edges;

% Water depth d is loaded from grid file above
d=depth;                   

% Initial condition
initial_condition;

% Set up pointers
setup_pointers(Ni,Nj);

% Allocate space for all variables
Rh = zeros(Ni,Nj);
Su = zeros(Ni+1,Nj);
Sv = zeros(Ni,Nj+1);
utilde = zeros(Ni+1,Nj);
vtilde = zeros(Ni,Nj+1);
Huwx = zeros(Ni+1,Nj);
Huwy = zeros(Ni,Nj+1);
Ap = ones(Ni+1,Nj);
Am = ones(Ni+1,Nj);
Bp = ones(Ni,Nj+1);
Bm = ones(Ni,Nj+1);
Ux = zeros(Ni+1,Nj);
Uy = zeros(Ni,Nj+1);
C2x = zeros(Ni+1,Nj);
C2y = zeros(Ni,Nj+1);
a = zeros(Ni+1,Nj);
b = zeros(Ni,Nj+1);
v_uface = zeros(Ni+1,Nj);
u_vface = zeros(Ni,Nj+1);

% Load in h data file if stored from previous run
if(LOAD_H)
  if(~exist('hstore','var'))
    load(hdatafile);
  end
end

t_start=tic;
% Main time loop
for n=1:nsteps
  u_old=u;
  v_old=v;
  h_old=h;

  % Get u on the v faces and v on the u faces
  [u_vface,v_uface]=ufaces(u,v);

  % Eulerian-Lagrangian advection
  if(ELM)
      %[u,v]=elm_p(xu,yu,xv,yv,xc,yc,d,H_small,u,v,u_vface,v_uface,dt,xpoly,ypoly);
      [u,v]=elm(xu,yu,xv,yv,xc,yc,d,H_small,u,v,u_vface,v_uface,dt);
    if(VERBOSE)
      fprintf('Completed ELM traceback integration.\n');
    end      
  elseif((exist('ADVECTION','var') & ADVECTION) | ...
         (exist('DIFFUSION') & DIFFUSION))
      if(~exist('Cu','var') | ~exist('Cv','var'))
          Cu=[];
          Cv=[];
      end
      [u,v,Cu,Cv]=advection_diffusion(u,v,dx,dy,Cu,Cv,n,dt,nu,ADVECTION,DIFFUSION);
  end

  % Total water depth is independent of h if the problem is linear
  if(LINEAR)
    H = d + hoffset;

    [Huwx,Huwy]=fluxfaceheights(H,d-zb,hoffset*ones(Ni,Nj),u,v);
  else
    H = d+h-zb; 
    
    % Correct for small depths (CWC problems?)
    H(H<H_small)=H_small;

    [Huwx,Huwy]=fluxfaceheights(H,d-zb,h,u,v);
  end
  
  [Hdragx,Hdragy]=dragfaceheights(H);

  % Provisional value of u for use in the drag term inside the domain
  utilde(inu_ij)=u(inu_ij)-g*dt/(2*dx)*(h(in_ij)-h(in_im1j));
  vtilde(inv_ij)=v(inv_ij)-g*dt/(2*dy)*(h(in_ij)-h(in_ijm1));

  % Provisional u on v faces and v on ufaces
  [utilde_vface,vtilde_uface]=ufaces(utilde,vtilde);

  % Magnitude of velocity on faces
  Ux=sqrt(utilde.^2+vtilde_uface.^2);
  Uy=sqrt(utilde_vface.^2+vtilde.^2); 

  % Set boundary velocities to zero
  [Ux,Uy]=faces_to_zero(Ux,Uy,outu_ij,outv_ij,outu_ip1j,outv_ijp1,markuface,markvface);

  % Drag coefficient is function of roughness and depth
  % Buffer drag and buffer height are needed to increase drag when
  % flow becomes very shallow.
  %  Cdmax=(1/kappa*(log(2)+0.5-1)).^(-2);
  %  H_buffer=2*z0;
  Cdmax=100;    % Increase drag when H<H_buffer to avoid
                % wetting/drying problems
  H_buffer=0.05;
  Cdx=(1/kappa*(log(Hdragx/z0)+z0./Hdragx-1)).^(-2);
  Cdy=(1/kappa*(log(Hdragy/z0)+z0./Hdragy-1)).^(-2);
  Cdx(Hdragx<H_buffer)=Cdmax;
  Cdy(Hdragy<H_buffer)=Cdmax;

  if(exist('DRAGPATCH','var') & DRAGPATCH)
      [Cdx,Cdy]=dragpatch(Cdx,Cdy,xu,yu,xv,yv);
  end

  % Frictional term within the domain
  Ap(inu_ij)=1+(2*Cdx(inu_ij).*dt.*theta.*Ux(inu_ij))./(H(in_im1j)+H(in_ij));
  Am(inu_ij)=1-(2*Cdx(inu_ij).*dt.*(1-theta).*Ux(inu_ij))./(H(in_im1j)+H(in_ij));
  Bp(inv_ij)=1+(2*Cdy(inv_ij).*dt.*theta.*Uy(inv_ij))./(H(in_ij)+H(in_ijm1));
  Bm(inv_ij)=1-(2*Cdy(inv_ij).*dt.*(1-theta).*Uy(inv_ij))./(H(in_ij)+H(in_ijm1));

  % Spatially-uniform wind stress
  [U10,taux,tauy]=windstress(n,dt);

  % Explicit side of momentum on u- and v-faces
  Su(inu_ij)=Am(inu_ij)./Ap(inu_ij).*u(inu_ij)...
      -g*(1-theta)*dt/dx./Ap(inu_ij).*(h(in_ij)-h(in_im1j))...
      +2*dt*taux./(Ap(inu_ij).*(H(in_ij)+H(in_im1j))) ...
      +f*v_uface(inu_ij).*dt./Ap(inu_ij);
  Sv(inv_ij)=Bm(inv_ij)./Bp(inv_ij).*v(inv_ij)...
      -g*(1-theta)*dt/dy./Bp(inv_ij).*(h(in_ij)-h(in_ijm1))...
      +2*dt*tauy./(Bp(inv_ij).*(H(in_ij)+H(in_ijm1)))...
      -f*u_vface(inv_ij).*dt./Bp(inv_ij);

  % Set the values of u, v, and h at the boundaries
  [u,v]=boundaryconditions_face([],[],[],u,v,Huwx*dy,Huwy*dx,xu,yu,xv,yv,n,boundaryregion_type2,Qb_type2,[],[],[],'Q');
  %[u,v]=boundaryconditions_face([],[],[],u,v,Huwx*dy,Huwy*dx,xu,yu,xv,yv,n,boundaryregion_type2,[],ub_type2,vb_type2,[],'u');
  h=boundaryconditions_cell(h,xc,yc,n,boundaryregion_type3,hb_type3);

  % Set the explicit side to 0 at boundary edges
  Su(markuface~=0)=0;
  Sv(markvface~=0)=0;
  Su(markuface==2)=u(markuface==2);
  Sv(markvface==2)=v(markvface==2);

  % Vectors for the pentadiagonal
  C2x = g*Huwx*dt^2/dx^2;
  C2y = g*Huwy*dt^2/dy^2;
  aih = theta^2*C2x./Ap;
  bjh = theta^2*C2y./Bp;
  
  % Explicit side of the free-surface eqution
  Rh = h - theta*dt/dx*(Huwx(2:end,:).*Su(2:end,:)-Huwx(1:end-1,:).*Su(1:end-1,:)) ...
      -(1-theta)*dt/dx*(Huwx(2:end,:).*u_old(2:end,:)-Huwx(1:end-1,:).*u_old(1:end-1,:)) ...
      - theta*dt/dy*(Huwy(:,2:end).*Sv(:,2:end)-Huwy(:,1:end-1).*Sv(:,1:end-1)) ...
      -(1-theta)*dt/dy*(Huwy(:,2:end).*v_old(:,2:end)-Huwy(:,1:end-1).*v_old(:,1:end-1));
  
  % set nonused cell to be 0
  Rh(cellmark==0)=0;
  
  % Set right-hand side of free-surface equation at boundary
  % cells to the boundary values (which are set in h)
  Rh(cellmark==3)=h(cellmark==3);
 
  h_old = h;
  % Compute the free-surface with the 2D conjugate-gradient
  % solver
  if(LOAD_H)
    h=hstore{n};
  else
    if(VERBOSE) 
      fprintf('CG Solve: ');
      tic;
    end
    
    [h,error] = cg2('fs2d',Rh,h,tol,nmax);

    if(VERBOSE) 
      fprintf('completed after %d iterations in %.3f s\n',length(error),toc);
    end
  end
  if(STORE_H)
    hstore{n}=h;
  end

  % Update velocities with new free-surface pressure gradient
  u(inu_ij)=Su(inu_ij)-g*theta*dt./Ap(inu_ij).*(h(in_ij)-h(in_im1j))/dx;
  v(inv_ij)=Sv(inv_ij)-g*theta*dt./Bp(inv_ij).*(h(in_ij)-h(in_ijm1))/dy;

  % Set boundary edges to zero (except type 2 edges)
  [u,v]=faces_to_zero(u,v,outu_ij,outv_ij,outu_ip1j,outv_ijp1,markuface,markvface);

  % Since the cg solve is not exact, need to update the free
  % surface using the velocity field to ensure that the continuity
  % equation is exact.  This ensures exact CWC regardless of the
  % tolerance of the cg solver.
  h = h_old - theta*dt/dx*(Huwx(2:end,:).*u(2:end,:)-Huwx(1:end-1,:).*u(1:end-1,:)) ...
      -(1-theta)*dt/dx*(Huwx(2:end,:).*u_old(2:end,:)-Huwx(1:end-1,:).*u_old(1:end-1,:)) ...
      - theta*dt/dy*(Huwy(:,2:end).*v(:,2:end)-Huwy(:,1:end-1).*v(:,1:end-1)) ...
      -(1-theta)*dt/dy*(Huwy(:,2:end).*v_old(:,2:end)-Huwy(:,1:end-1).*v_old(:,1:end-1));
  h(cellmark==3)=h_old(cellmark==3);

  % Set velocity on type 3 edges to interior values - this is only
  % needed when advection or diffusion of momentum are used
  if((exist('ADVECTION','var') & ADVECTION) | ...
         (exist('DIFFUSION') & DIFFUSION))
      u(u_type3u)=u(u_type3u_inner1);
      v(v_type3u)=v(v_type3u_inner1);
      u(u_type3v)=u(u_type3v_inner1);
      v(v_type3v)=v(v_type3v_inner1);
      %u(u_type3u)=2*u(u_type3u_inner1)-u(u_type3u_inner2);
      %v(v_type3u)=2*v(v_type3u_inner1)-v(v_type3u_inner2);
      %      u(u_type3v)=2*u(u_type3v_inner1)-u(u_type3v_inner2);
      %      v(v_type3v)=2*v(v_type3v_inner1)-v(v_type3v_inner2);
  end

  % Erosion/Deposition into rhs_C which is used in SSC and bedload functions
  if(((exist('EXNER','var') & EXNER) | ...
      (exist('SSC','var') & SSC)))
      rhs_C=zeros(Ni,Nj);
      
      if(exist('EROSION','var') & EROSION)
          rhs_C=rhs_C+erosion(n,dt,C,u_old,v_old,zb,xc,yc,h_old,d,cellmark);
      end
      if(exist('DEPOSITION','var') & DEPOSITION)
          rhs_C=rhs_C-ws*C;
      end
  end

  % Bedload
  if((exist('BEDLOAD','var') & BEDLOAD) & (exist('HELICAL','var') & HELICAL))
      if(exist('HELICAL_ADVECTION','var') & HELICAL_ADVECTION)
          [tan_deltas_0,rhs_helical] = helical_model_rhs(u,v,dx,dy,H,tan_deltas);
          tan_deltas=scalartransport(u_old,v_old,u,v,xc,yc,xu,yu,xv,yv,tan_deltas,...
                                     Huwx,Huwy,dx,dy,dt,n,theta,h,d,h_old+d,L,W,rhs_helical,...
                                     boundaryregion_type2,tan_deltas_b_type2,...
                                     boundaryregion_type3, ...
                                     tan_deltas_b_type3);
      else
          tan_deltas=tan_deltas_0;
      end
  else
      tan_deltas=zeros(Ni,Nj);
  end
  if(exist('EXNER','var') & EXNER)
      zb = exner(u_old,v_old,zb,dx,dy,dt,H,tan_deltas,C,-rhs_C,n);
  end
  if(exist('SSC','var') & SSC)
      C=scalartransport(u_old,v_old,u,v,xc,yc,xu,yu,xv,yv,C,...
                        Huwx,Huwy,dx,dy,dt,n,theta,h,d,h_old+d,L,W,rhs_C,...
                        boundaryregion_type2,Cb_type2,...
                        boundaryregion_type3,Cb_type3);

      % Check for CWC here
      %      fprintf('min(C-C0) = %.3e, max(C-C0) =
      %      %.3e\n',min(min(C-100)),max(max(C-100)));
  end

  % Output data at specified points given by stationlocations.x and stationlocations.y
  if(~isempty(vars_to_output))
      outputswe;
  end

  % output free-surface results
  if(PLOT & (rem(n,ntout)==0 | n==1))
      plotswe(xc,yc,u,v,h,C,zb,d,cellmark,n,nsteps,dt,vars_to_plot,MOVIE);
  end

  % Report Cmax and time remaining  
  if(~rem(n,ntprog) | n==1 | n==nsteps)
      progress(u,v,C,t_start,n,dt,dx,dy,nsteps,Cmax_allowed,C0_max_allowed);
  end

  % User-defined function for processing data on the fly
  if(exist('process_data.m','file'))
      process_data;
  end
end 

if(STORE_H)
  save(hdatafile,'hstore');
end







