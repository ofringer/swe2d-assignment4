%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% CG2 2-d Conjugate gradient solver
%   [x,error] = CG2('L',b,tol,nmax) solves the 2-d linear
%   system defined by Lx=b with the conjugate gradient
%   algorithm, where the function L.m takes
%   as its argument a column vector the length of b and
%   returns L(x).  The operator must be symmetric and
%   positive definite.  The tolerance is specified by
%   tol and nmax is the maximum number of iterations, and
%   the error history is returned in the error vector.
%
% Oliver Fringer
% Stanford University
% 18 May 2004
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,error] = cg2(L,b,x0,tol,nmax)

  global in_ij

  error = 0;
  x = x0;
  r = b-eval([L,'(x)']);
  p = r;

  if(norm(b,'fro')==0)
    x = b;
    error = 0;
    return;
  end
  
  n=2;
  bnorm=norm(b(:));
  rnorm=norm(r(:));
  error(1)=inf;
  while(error(n-1)>tol & n<nmax)

    error(n)=rnorm/bnorm;
    beta = rnorm^2;
    z = eval([L,'(p)']);
    alpha = beta/(p(in_ij)'*z(in_ij));
    x = x + alpha*p;
    r = r - alpha*z;
    rnorm = norm(r(in_ij));
    beta = rnorm^2/beta;
    p = r + beta*p;

    n = n+1;
  end

  if(n>=nmax)
    fprintf('Warning!  Iteration not converging after %d iterations!\n',n);
  end
