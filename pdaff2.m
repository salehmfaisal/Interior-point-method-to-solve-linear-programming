%
%  [x,y,z]=pdaff(A,b,c,x0,y0,z0,tau,maxiter)
%
%  Inputs:
%     A,b,c        The problem, min c'*x, Ax=b, x>=0
%     x0,y0,z0     Initial solution, with x0,z0>0
%
%  Outputs:
%     x,y,z        An optimal solution.
%
function [x,y,z]=pdaff2(A,b,c,x0,y0,z0,tau,maxiter)
%
% Get the problem size.
%
[m,n]=size(A);
%
% First, set tau if it wasn't given.
%
if (nargin < 7)
  tau=1.0e-6;
end
%
% Also set maxiter if it wasn't given.
%
if (nargin < 8)
  maxiter=100;
end
%
% Initialize the loop 
%
x=x0;
y=y0;
z=z0;
ret=0;
iter=0;
%
% Main loop, do primal-dual affine method iterations.
%
while (iter<maxiter)
  if (((norm(A*x-b)/(1+norm(b))) < tau) & ...
      ((norm(A'*y+z-c)/(1+norm(c)) < tau) & ...
      (abs(c'*x-b'*y)/(1+abs(c'*x))) < tau))
    return
  else
%
%  Compute the step.
%
    delta=[zeros(n,n) A' eye(n); A zeros(m,m) zeros(m,n); ...
           diag(z) zeros(n,m) diag(x)]\[-A'*y-z+c; -A*x+b; -x.*z];
    deltax=delta(1:n);
    deltay=delta(n+1:n+m);
    deltaz=delta(n+m+1:2*n+m);
%
%  Find the maximum alphap
%
    alphap=1;
    for i=1:n,
      if (x(i)+alphap*deltax(i) < 0),
        alphap=-x(i)/deltax(i);
      end;
    end;
    alphap=alphap*.99;
%
% Find the maximum alphad
%
    alphad=1;
    for i=1:n,
      if (z(i)+alphad*deltaz(i) < 0),
        alphad=-z(i)/deltaz(i);
      end;
    end;
    alphad=alphad*.99;
%
% Take the step.
%
    x=x+alphap*deltax;
    y=y+alphad*deltay;
    z=z+alphad*deltaz;
%
% Output some useful information
%
    iter=iter+1;
%    alphap
%    alphad
%    pinfeas=norm(A*x-b)/(1+norm(b))
%    dinfeas=norm(A'*y+z-c)/(1+norm(c))
%    pobj=c'*x
%    dobj=b'*y
     x(1:2)
%    pause;

  end;
end;
