function rhs_C=erosion(n,dt,C,u,v,zb,xc,yc,h,d,cellmark)

% Fetch and depth along the fetch are global since they are only
% computed on the first time step
global F Df

getvariables;

% Grid size
[Ni,Nj] = size(d);

% Water depth
H = h+d;
H0 = hoffset+d;

if(n==1 & FETCH)
    [F,Df] = fetch(thetaW,xc,yc,H,cellmark);
end 

% Wind stress
[U10,taux_wind,tauy_wind]=windstress(n,dt);

% Example use of dispersion function with random frequencies
omega_wave = rand(Ni,Nj);
    
k = dispersion(omega_wave,g,H);      % Wavenumber from linear dispersion relation

% Need to set rhs_C to erosion
rhs_C = zeros(Ni,Nj);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Dispersion relation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function k = dispersion(omega,g,H)

Y = omega.^2.*H/g;
X = zeros(size(Y));
X(Y>=2) = Y(Y>=2); 

x0 = linspace(0,3,100);
y0 = x0.*tanh(x0);

Y(Y<=0) = 0;

if(min(min(Y))<min(y0) | max(max(Y(Y<2)))>max(y0))
  fprintf('Max y0 is %f\n',max(y0));		       
  error('Y out of bounds');
end
if(myisnan(Y))
  error('Y is nan');
end
X(Y<2) = interp1(y0,x0,Y(Y<2));
if(myisnan(X))
  error('X is nan');
end
k = X./H;
k(Y<=0) = 0;
k(H<=0) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% My isnan
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = myisnan(x)
 
  if(~isempty(find(isnan(x))))
    y=true;
  else
    y=false;
  end
    