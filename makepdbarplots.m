%
% Load the second example problem.  
%
load ex2.mat
x=x0;
y=y0;
z=z0;
epsilon=1.0e-8;
sols=zeros(2,8);
sols(:,1)=x(1:2);
for i=1:7
  [x,y,z]=simplepd(A,b,c,x,y,z,epsilon,1);
  sols(:,i+1)=x(1:2);
end
%
% Plot the solutions.
%
figure(1);
clf;
%
% First, plot the feasible region.
%

hold on
plot([0 0],[0 5],'k-');
plot([0 5],[0 0],'k-');
grid
%
% Plot the constraints as solid lines.
%
xs=0:5;
ys=4*ones(size(xs));
plot(xs,ys,'k-');
ys=5*ones(size(xs))-xs;
plot(xs,ys,'k-');
ys=7*ones(size(xs))-2*xs;
plot(xs,ys,'k-');
ys=9.5*ones(size(xs))-3*xs;
plot(xs,ys,'k-');
ys=12.3333333*ones(size(xs))-4*xs;
plot(xs,ys,'k-');
%
% Now, plot the solutions.
% 
plot(sols(1,:),sols(2,:),'ko--');
axis([0 5 0 5]); 
print -deps pdbarfig1.eps
%eps2pdf('pdbarfig1.eps');
%
% Load the example problem and do it again, starting at a different
% point.
%
load ex2.mat
x=x0;
x(1)=0.1;
x(2)=4;
y=y0;
z=z0;
epsilon=1.0e-8;
sols=zeros(2,10);
sols(:,1)=x(1:2);
for i=1:9
  [x,y,z]=simplepd(A,b,c,x,y,z,epsilon,1);
  sols(:,i+1)=x(1:2);
end
%
% Plot the solutions.
%
figure(2);
clf;
%
% First, plot the feasible region.
%
hold on
plot([0 0],[0 5],'k-');
plot([0 5],[0 0],'k-');
grid
%
% Plot the constraints as solid lines.
%
xs=0:5;
ys=4*ones(size(xs));
plot(xs,ys,'k-');
ys=5*ones(size(xs))-xs;
plot(xs,ys,'k-');
ys=7*ones(size(xs))-2*xs;
plot(xs,ys,'k-');
ys=9.5*ones(size(xs))-3*xs;
plot(xs,ys,'k-');
ys=12.3333333*ones(size(xs))-4*xs;
plot(xs,ys,'k-');
%
% Now, plot the solutions.
% 
plot(sols(1,:),sols(2,:),'ko--');
axis([0 5 0 5]); 
print -deps pdbarfig2.eps
%eps2pdf('pdbarfig2.eps');
