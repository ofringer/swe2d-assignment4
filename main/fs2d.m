%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FS2DNEW Return the left-hand side operator for the conjugate 
%      gradient solution of the shallow water equations.  This 
%      operator is given by
% 
%      L(h_{i,j}) = -a_{i,j}h_{i-1,j} - b_{i,j}h_{i,j-1}
%          +(1 + a_{i,j} + a_{i+1,j} + b_{i,j} + b_{i,j+1})h_{i,j}
%                   -a_{i+1,j}h_{i+1,j} - b_{i,j+1}h_{i,j+1}
%      
%      except in this m-file we use biph,bimh,aiph,aimh to represent
%      the values of the coefficients at the faces.
%      Here the new program will work with swe2dnew to return left-hand
%      side accoring to cell location and bc type
%      
% Oliver Fringer
% Yun Zhang
% Stanford University
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = fs2d(x)
  
  global aih bjh in_ij inu_ij in_im1j inu_ip1j in_ip1j inv_ij in_ijm1 ...
      inv_ijp1 in_ijp1 cellmark outu_ij outv_ij outu_ip1j outv_ijp1
 
  [Ni,Nj]=size(x);
  y=x;

  aih(Ni+1,:)=0;
  bjh(:,Nj+1)=0;

  aih(outu_ij)=0;
  bjh(outv_ij)=0;

  aih(outu_ip1j)=0;
  bjh(outv_ijp1)=0;

  y(in_ij) = x(in_ij) + aih(inu_ij).*(x(in_ij)-x(in_im1j)) - aih(inu_ip1j).*(x(in_ip1j)-x(in_ij))...
      +bjh(inv_ij).*(x(in_ij)-x(in_ijm1)) - bjh(inv_ijp1).*(x(in_ijp1)-x(in_ij));

  y(cellmark==3)=x(cellmark==3);
  y(cellmark==0)=0;
