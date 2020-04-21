x=[0 0 0 0]';
A=[5 2 1 1;2 6 2 1;1 2 7 2;1 1 2 8];
b=[29;31;26;19];
tol=10^-6;
d=diag(A);
M=diag(d);
r=b-A*x;
p=inv(M)*r;
k=0;
while (norm(b-A*x)>tol)
    q=A*p;
    a=(p'*r)/(p'*q);
    x=x+a*p
    r=r-a*q;
    p=inv(M)*r;
    k=k+1
end

