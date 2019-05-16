% switches for functions
LINEAR=false;           % Linear free surface, i.e. H = d, not H=h+d;
ELM=true;               % Whether or not to use ELM
VERBOSE=false;          % Verbose output
TRANSPORT=true;         % Compute scalar transport
TIDES=true;             % Whether or not to include tides

STORE_H=false;          % Store h for offline use
LOAD_H=false;           % Load h so don't need to use slow CG solve
hdatafile='data/hdata.mat';  % Data where h is stored

% Gravity 
g = 9.81;

% Maximum allowable Courant number before exiting
Cmax_allowed = 2;

% Coriolis parameter
Tday = 86400;           % Seconds per day
Omega_day = 2*pi/Tday;  % diurnal frequency
phi = 37;               % Latitude in degreess
f=2*Omega_day*sin(phi*pi/180); % in radians

% Wind stress in x and y directions
tau_x0 = 0.0;
tau_y0 = 0.0;

% For CG solver
tol = 1e-5;           % Tolerance for conjugate gradient solver
nmax = 400;            % Maximum number of iterations for conjugate gradient solver

% Bottom roughness
z0 = 1e-4;
kappa = 0.41;

% Frequency of wind-stress forcing (in rad/s)
omega_0 = 0;

% For plotting (1000 m = 1 km);
km=1000; 

% For tidal constituents
tidefile='data/tidedata_point_reyes_2005';

% Time step in s
dt=0.05; % Curved Ch

tmax = 40;             % total simulation time
nsteps = fix(tmax/dt); % total number of time steps
tout = dt*10;          % how often to plot output (in hr)
tmax = nsteps*dt;      % reset tmax
theta = 0.55;          % implicitness parameter

% Mean sea-level is offset by 4 m relative to bathymetry datum
hoffset=0;

% Bathymetry is in depth
Dmin = -inf;         % Minimum allowable bottom depth (relative to
                     % MSL) to prevent wetting/drying (otherwise
                     % the calculation takes a lot longer).

% Smallest water depth (h+d) allowed before correction
H_small = 0.01;

% Plotting
PLOT=true;              % Plot output
vars_to_plot='qC';

% Boundary conditions
boundaryregion_type2{1} = [];[4 7 -1 1];
Qb_type2{1} = 0.053;
Cb_type2{1} = 1;
ub_type2{1} = 0;
vb_type2{1} = 0.7;

boundaryregion_type3{1} = [-1 3 -1 1];
Cb_type3{1} = 0;
hb_type3{1} = 0;

boundaryregion_type3{2} = [4 7 -1 1];
Cb_type3{2} = 1;
hb_type3{2} = 0.05;

