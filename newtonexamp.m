%
% First, plot the test function.
%
figure(1);
clf;
x=-2:0.01:2;
y=testfun(x);
plot(x,y,'k');
hold on
%
% Plot x, y axes and a grid.
%
grid on
plot([-2 2],[0 0],'k');
plot([0 0],[-6 12],'k');
xlabel('x');
ylabel('f(x)');
%
% Plot the tangent line at x=1.5.
%
x0=1.5;
fx0=testfun(x0);
%
% Plot the point.
%
plot(x0,0,'k+');
plot(x0,fx0,'k+');
%
% Get the tangent line and plot it.
%
fpx0=testfunp(x0);
ln=fx0+fpx0*(x-x0);
plot(x,ln,'k--');
axis([-2 2 -6 8]);
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
fx1=testfun(x1);
fpx1=testfunp(x1);
ln2=fx1+fpx1*(x-x1);
plot(x,ln2,'k--');
axis([-2 2 -6 8]);
plot(x1,fx1,'k+');
x2=x1-fx1/fpx1;
plot(x2,0,'k+');
%
% Add labels for x0, x1, x2.
%
text(x0,-0.7,'x_0');
text(x1,-0.7,'x_1');
text(x2,-0.7,'x_2');
%
% Print it out.
%
print -deps fig1.eps
