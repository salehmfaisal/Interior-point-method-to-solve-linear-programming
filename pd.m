function [x,y,z,iter]=pd(A,b,c,x0,y0,z0,tol,maxiter)
%
%  [x,y,z,iter]=pd(A,b,c,x0,y0,z0,tol,maxiter)
%
%  Implements a simple primal dual method with sigmak=1/2, using
%  sparse Cholesky factorization.
%
%  Inputs:
%     A,b,c        The problem, min c'*x, Ax=b, x>=0
%     x0,y0,z0     Initial solution, with x0,z0>0
%     tol          Convergence tolerance (default 1.0e-6.)
%     maxiter      Maximum number of iterations to run (default 100.)
%
%  Outputs:
%     x,y,z        An optimal solution.
%     iter         Iteration count.
%

%
% Get the problem size.
%
[m,n]=size(A);
%
% First, set tol if it wasn't given.
%
if (nargin < 7)
  tol=1.0e-6;
end
%
% Also set maxiter if it wasn't given.
%
if (nargin < 8)
  maxiter=100;
end

%
% Compute the best permutation of A*Zi*X*A'.  Note that this only needs
% to be done once, since the permutation will never change.
% Furthermore, this is indepedent of the actual values of Zi and X, 
% so we can simply use A*I*I*A'=A*A'.
%
M=A*A';
p=symamd(M);
%
% Initialize the loop 
%
x=x0;
y=y0;
z=z0;
deltay=zeros(size(y));
deltax=zeros(size(x));
deltaz=zeros(size(z));
ret=0;
iter=0;
while (iter<maxiter)
  if (((norm(A*x-b)/(1+norm(b))) < tol) & ...
      ((norm(A'*y+z-c)/(1+norm(c)) < tol) & ...
      (abs(c'*x-b'*y)/(1+abs(c'*x))) < tol))
    iter
    return
  else
%
% Compute mu.
%
  mu=x'*z/n;
%
% Pick sigma=1/2.
%
  sigma=1/2;
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
% Compute the system matrix and rhs.
%
  M=A*Zi*X*A';
  rhs=-M*y+A*X*Zi*c-sigma*mu*A*Zi*e+b-A*x;
  
%
% Uncomment the following line to see the condition number of M.
% As we approach an optimal solution, M typically becomes badly conditioned.
%
%  fprintf(1,'cond(M)=%f\n',cond(full(M)))
%
%  Apply the permutation.
%
  M=M(p,p);
  rhs=rhs(p);
%
% Compute the deltay step.
%
   R=chol(M);
   v=R'\rhs;
   deltay(p)=R\v;
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
% Output some useful information.
%
    iter=iter+1
    alphap;
    alphad;
    pinfeas=norm(A*x-b)/(1+norm(b));
    dinfeas=norm(A'*y+z-c)/(1+norm(c));
    pobj=c'*x
    dobj=b'*y
  end;
end;
