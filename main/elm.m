%
% Compute the velocity vectors at the location of the Lagrangian
% tracebacks.
%
% Oliver Fringer
% Yun Zhang
% Stanford University
%
function [uminus,vminus]=elm(xu,yu,xv,yv,xc,yc,d,H_small,u,v,u_vface,v_uface,dt)

global_pointers;

% Calculate the Lagrangian traceback from the u-faces
xminus_uface = xu - dt*u;
yminus_uface = yu - dt*v_uface;

% Determine intersection of tracebacks with boundaries
% Need to create xbox and ybox to create boundary polygon
%[xi, yi] = polyxpoly(x, y, xbox, ybox, 'unique');

% Only update the velocities for which the traceback is
% within the domain, or where the depth is greater than H_small
in = incheck(xminus_uface,yminus_uface,xc',yc',d',H_small);
xminus_uface(~in) = xu(~in);
yminus_uface(~in) = yu(~in);
uminus = interp2(xu',yu',u',xminus_uface,yminus_uface);
uminus(isnan(uminus))=0;

% Calculate the Lagrangian traceback from the v-faces
xminus_vface = xv - dt*u_vface;
yminus_vface = yv - dt*v;

% Only update the velocities for which the traceback is
% within the domain, or where the depth is greater than H_small
in = incheck(xminus_vface',yminus_vface',xc',yc',d',H_small);
xminus_vface(~in) = xv(~in);
yminus_vface(~in) = yv(~in);
vminus = interp2(xv',yv',v',xminus_vface,yminus_vface);
vminus(isnan(vminus))=0;

% Don't do ELM on boundary edges.
uminus(markuface>0)=u(markuface>0);
vminus(markvface>0)=v(markvface>0);

