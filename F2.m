function y=F2(x)
%
% To avoid getting a complex value, we'll just return NaN if x<=0.
%
if (x>0)
  y=log(x)-1;
else
  y=NaN;
end