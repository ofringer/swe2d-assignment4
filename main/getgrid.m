function [dx,dy,Ni,Nj,depth,markuface,markvface,cellmark]=getgrid(gridfile)


%[dx,dy,Ni,Nj,depth,markuface,markvface,cellmark]=getrectangulargrid();

%load(gridfile);

[dx,dy,Ni,Nj,depth,markuface,markvface,cellmark]=curvedchannel();
