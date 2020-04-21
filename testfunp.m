%
% yp=testfunp(x)
%
% Derivative of a simple test function for Newton's method.  We compute 
%
%   y=x^3+x^2-1
%   yp=3x^2+2x
%
function yp=testfunp(x)
yp=3*x.^2+2*x;