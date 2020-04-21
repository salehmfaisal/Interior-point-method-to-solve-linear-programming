%
% L=mycholinf(M)
%
% Computes the Cholesky-Infinity factorization of a symmetric and
% positive definite matrix $M$.  
%
function L=mycholinf(M)
n=size(M,1);
L=zeros(n,n);
for j=1:n
  L(j,j)=M(j,j);
  for k=1:(j-1)
    L(j,j)=L(j,j)-L(j,k)^2;
  end
  if (L(j,j)>0)
    L(j,j)=sqrt(L(j,j));
  else
    L(j,j)=+Inf;
  end
  for i=(j+1):n
    L(i,j)=M(i,j);
    for k=1:(j-1)
      L(i,j)=L(i,j)-L(i,k)*L(j,k);
    end
    L(i,j)=L(i,j)/L(j,j);
  end
end 