% SSC is used for scalar transport when there is no settling or
% erosion
SSC=true;
EROSION=true;
DEPOSITION=true;
EXNER=true;
FETCH=true;

% Whether or not to use TIDES
TIDES=false;
FLUXES=true;

% Linear free surface, i.e. H = d, not H=h+d;
LINEAR=false;                             

% ELM is not needed for this assignment
ELM=true;

% Verbose output
VERBOSE=false;                            

% Store h for offline use
STORE_H=false;

% Load h so don't need to use slow CG solve
LOAD_H=false;

% Whether or not to create a movie
MOVIE=false;

% Data where h is stored
hdatafile='data/hdata.mat';               

% Gravity 
g = 9.81;                                 

% Sediment parameters
ws = 5e-4;

% Seconds per day
Tday = 86400;                             

% Coriolis parameter
Omega_day = 2*pi/Tday;                    
phi = 37;                                 
f=2*Omega_day*sin(phi*pi/180);            

% Wind parameters
thetaW = -25*pi/180; % 90 or -25
Uw0 = 4;
Cdw = 0.0005;
rho_air = 1.23;
rho0 = 1000;

% Frequency of wind-stress forcing (in rad/s)
omega_0 = 0;

% Kinematic viscosity
nu = 1e-6;

% Tolerance for conjugate gradient solver
tol = 1e-4;                               

% Maximum number of iterations for conjugate gradient solver
nmax = 400;                               

% Bottom roughness
z0 = 1e-4;                                

% Von karman constant
kappa = 0.41;                             

% For plotting (1000 m = 1 km);
km=1000;                                  

% Time step size
dt=360;

% Maximum allowable Courant number before exiting
Cmax_allowed = 1;                         

% Maximum allowable SSC (to detect blowups for scalar transport)
C0_max_allowed=10000;

% Total simulation time
max_time = 8*Tday; 

% total number of time steps
nsteps = fix(max_time/dt);                    
tmax = nsteps*dt; % Do not change                         

% how often to plot output
ntout = fix(6*3600/dt);                     

% implicitness parameter
theta = 0.55;                             

% Mean sea-level is offset relative to bathymetry datum
hoffset=-4;                               

% How often to report progress
ntprog = 100;

% Minimum allowable bottom depth (relative to MSL) to prevent wetting/drying (otherwise
% the calculation takes a lot longer)
Dmin = -inf;                              
                                          
% Smallest water depth (h+d) allowed before correction
H_small = 0.1;                            

% Plot output
PLOT=true;

% List of variables to plot, C=sediment/scalar, u/v = velocity components, h=free surface, q=quiver plot:
vars_to_plot='q';

% Variables needed to output time series and the locations at which
% those time series are desired.
nout=3;
vars_to_output='uvChz';
x_transect=23000;
y_transect=15000;

% Indices to edges to compute fluxes - used in fluxes.m
flux_file = 'data/dumbarton_flux_indices.mat';
          
