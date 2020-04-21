%
% yp=testfun2p(x)
%
% Derivative of a simple test function for Newton's method.  We compute 
%
%   
function yp=testfun2p(x)
yp=exp(-x).*(-1+1./x+1./x.^2);
