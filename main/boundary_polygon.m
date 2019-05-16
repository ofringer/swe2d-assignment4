function [xb,yb]=boundary_polygon(xu,yu,xv,yv,markuface,markvface)

xb = xu(markuface~=0);

x(mark==0)=[];
y(mark==0)=[];
mark(mark==0)=[];

figure(1)
plot(x,y,'k.');
axis image;

xb=x;
yb=y;