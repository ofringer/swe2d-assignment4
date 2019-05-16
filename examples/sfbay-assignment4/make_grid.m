function [dx,dy,Ni,Nj,depth,markuface,markvface,cellmark]=make_grid()

load('data/southbay_grid');

% Bathymetry is in depth
hoffset=-4;
H = depth+hoffset;
Dmin = 1.0;               % Minimum allowable bottom depth (relative to
                        % MSL) to prevent wetting/drying (otherwise
                        % the calculation takes a lot longer).
H(H<Dmin & depth~=0)=Dmin;
depth = H-hoffset;

