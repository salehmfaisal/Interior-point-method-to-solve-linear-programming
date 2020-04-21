%
%  [x,y,z]=pdhalf(A,b,c,x0,y0,z0)
%
%  Implements a simple primal dual method with sigmak=1/2.
%
%  Inputs:
%     A,b,c        The problem, min c'*x, Ax=b, x>=0
%     x0,y0,z0     Initial solution, with x0,z0>0
%
%  Outputs:
%     x,y,z        An optimal solution.
%
function [x,y,z]=pdhalf(A,b,c,x0,y0,z0)
x=x0;
y=y0;
z=z0;
[m,n]=size(A);
ret=0;
iter=0;
alphap=1;
alphad=1;
while (ret == 0)
  if (((norm(A*x-b)/(1+norm(x))) < 1.0e-6) & ...
      ((norm(A'*y+z-c)/(1+norm(y)+norm(z))) < 1.0e-6) & ...
      (abs(c'*x-b'*y)/(1+abs(c'*x)) < 1.0e-6))
    ret=1;
  else
%
% Compute mu.
%
  mu=x'*z/n;
%
% Pick sigma=1/2 or sigma=1.  Use sigma=1 if the last step was short.
%
  sigma=1/n;

  if (alphap+alphad < .5)
    sigma=1;
  end;

%
%  Compute the step.
%
    delta=[zeros(n,n) A' eye(n); A zeros(m,m) zeros(m,n); ...
           diag(z) zeros(n,m) diag(x)]\[-A'*y-z+c; -A*x+b; -x.*z + sigma*mu*ones(n,1)];
    deltax=delta(1:n);
    deltay=delta(n+1:n+m);
    deltaz=delta(n+m+1:2*n+m);
%
%  Find the maximum alphap
%
    alphap=1/.99;
    for i=1:n,
      if (x(i)+alphap*deltax(i) < 0),
        alphap=-x(i)/deltax(i);
      end;
    end;
    alphap=alphap*.99;
%
% Find the maximum alphad
%
    alphad=1/.99;
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
    iter=iter+1
    alphap
    alphad
    pinfeas=norm(A*x-b)/(1+norm(x))
    dinfeas=norm(A'*y+z-c)/(1+norm(y)+norm(z))
    pobj=c'*x
    dobj=b'*y

    pause;
    
  end;
end;
