%
% x=newton(fname,Jname,x0,maxiter,tol)
%
function x=newton(fname,Jname,x0,maxiter,tol)
x=x0;
F=feval(fname,x);
J=feval(Jname,x);
iter=0;
while ((iter < maxiter) & (norm(F)>tol)),
    deltax=(J)\(-F);
    x=x+deltax
    F=feval(fname,x);
    J=feval(Jname,x);
    iter=iter+1;
end;
