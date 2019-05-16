%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Fetch function
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [F,Df] = fetch(thetaW,xc,yc,d,cellmark)

  d(cellmark==0)=0;
  L = max(max(xc))-min(min(xc));
  W = max(max(yc))-min(min(yc));
  [is,js]=find(cellmark~=0);
  Df = zeros(size(cellmark));
  F = zeros(size(cellmark));

  F = zeros(size(xc));
  for m=1:length(is)
    xs = xc(is(m),js(m));
    ys = yc(is(m),js(m));
    xe = xs-L*cos(thetaW);
    ye = ys-L*sin(thetaW);
    xl = linspace(xs,xe,100);
    yl = linspace(ys,ye,100);
    rl = sqrt((xl-xs).^2+(yl-ys).^2);
    dl = interp2(xc',yc',d',xl,yl);
    dl(isnan(dl))=0;

    imin=find(dl==0);
    F(is(m),js(m)) = rl(min(imin));
    Df(is(m),js(m)) = mean(dl(1:imin));
  end  
