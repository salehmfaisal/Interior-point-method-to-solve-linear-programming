%
% First, plot the test function.
%
figure(2);
clf;
x=0:0.01:20;
y=testfun2(x);
plot(x,y,'k');
hold on
%
% Plot x, y axes and a grid.
%
grid on
plot([0 15],[0 0],'k');
plot([0 0],[-10 5],'k');
xlabel('x');
ylabel('f(x)');
%
% Plot the tangent line at x0=2.
%
x0=2;
fx0=testfun2(x0);
%
% Plot the point.
%
plot(x0,0,'k+');
plot(x0,fx0,'k+');
%
% Get the tangent line and plot it.
%
fpx0=testfun2p(x0);
ln=fx0+fpx0*(x-x0);
plot(x,ln,'k--');
%
% Plot the next point.
%
x1=x0-fx0/fpx0;
plot(x1,0,'k+');
%
% and, do the second iteration.
%
%
% Get the tangent line and plot it.
%
fx1=testfun2(x1);
fpx1=testfun2p(x1);
ln2=fx1+fpx1*(x-x1);
plot(x,ln2,'k--');
axis([0 7 -0.5 0.5]);
plot(x1,fx1,'k+');
x2=x1-fx1/fpx1;
plot(x2,0,'k+');
%
% Add labels for x0, x1, x2.
%
text(x0,-0.05,'x_0');
text(x1,-0.05,'x_1');
text(x2,-0.05,'x_2');
%
% Print it out.
%
print -deps fig2.eps
