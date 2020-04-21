%
%  [x,y,z]=pdhalf2(A,b,c,x0,y0,z0)
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
%  This version of the code is designed to use the Cholesky factorization
%  and exploit sparsity in the A matrix.  
%
function [x,y,z]=pdhalf(A,b,c,x0,y0,z0)
%
% Initialize some variables.
%
x=x0;
y=y0;
z=z0;
[m,n]=size(A);
ret=0;
iter=0;
alphap=1;
alphad=1;
deltay=zeros(m,1);
deltaz=zeros(n,1);
deltax=zeros(n,1);
%
%  The main loop.
%
while (ret == 0)
%
% First, check for convergence.
%
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


%
%  X=diag(x) and Z=diag(z).  While we're here, compute inv(X) and inv(Z).
%  Also compute a vector of all ones.
%
  X=spdiags(x,0,n,n);
  Z=spdiags(z,0,n,n);
  Xi=inv(X);
  Zi=inv(Z);
  e=ones(n,1);
%
% Compute the best permutation of A*Zi*X*A'.  Note that this only needs
% to be done once, since the permutation will never change. 
%
  if (iter == 0),
    P=symamd(A*Zi*X*A');
  end;
%
% Compute the system matrix and rhs.
%
  M=A*Zi*X*A';
  rhs=-M*y+A*X*Zi*c-sigma*mu*A*Zi*e+b-A*x;
%
%  Apply the permutation.
%
  M=M(P,P);
  rhs=rhs(P);
%
% Compute the deltay step.
%
   R=chol(M);
   foo=R'\rhs;
   deltay(P)=R\foo;
%
% Now, compute the deltax and deltaz steps.  
%
  deltax=Zi*X*(A'*(deltay+y)-c)+sigma*mu*Zi*e;
  deltaz=-Xi*Z*deltax-Z*e+sigma*mu*Xi*e;    

%
%  Find the maximum alphap
%
    alphap=1/0.99;
    for i=1:n,
      if (x(i)+alphap*deltax(i) < 0),
        alphap=-x(i)/deltax(i);
      end;
    end;
    alphap=alphap*.99;
%
% Find the maximum alphad
%
    alphad=1/0.99;
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
    alphap;
    alphad;
    pinfeas=norm(A*x-b)/(1+norm(x));
    dinfeas=norm(A'*y+z-c)/(1+norm(y)+norm(z));
    pobj=c'*x
    dobj=b'*y

    if (mod(iter,100) == 0),    
      pause;
    end;
  end;
end;
