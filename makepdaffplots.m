%
%
%
load ex1.mat
x=x0;
y=y0;
z=z0;
tau=1.0e-8;
sols=zeros(2,6);
sols(:,1)=x(1:2);
for i=1:5
  [x,y,z]=pdaff(A,b,c,x,y,z,tau,1);
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
plot([0 0],[0 4],'k-');
plot([0 4],[0 0],'k-');
grid
%
% Plot the constraints as solid lines.
%
xs=0:4;
ys=2*ones(size(xs));
plot(xs,ys,'k-');
ys=3*ones(size(xs))-xs;
plot(xs,ys,'k-');
plot([2 2],[0 4],'k-');
%
% Now, plot the solutions.
% 
plot(sols(1,:),sols(2,:),'ko--');
axis([0 4 0 4]); 
print -deps pdafffig1.eps
%eps2pdf('pdafffig1.eps');
%
% Load the second example problem.  
%
load ex2.mat
x=x0;
y=y0;
z=z0;
tau=1.0e-8;
sols=zeros(2,8);
sols(:,1)=x(1:2);
for i=1:7
  [x,y,z]=pdaff2(A,b,c,x,y,z,tau,1);
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
print -deps pdafffig2.eps
%eps2pdf('pdafffig2.eps');
%
% Load the example problem and do it again, starting at a different
% point.
%
load ex2.mat
x=x0;
x(1)=0.1
x(2)=4-1.0e-6;
y=y0;
z=z0;
tau=1.0e-8;
sols=zeros(2,20);
sols(:,1)=x(1:2);
for i=1:20
  [x,y,z]=pdaff2(A,b,c,x,y,z,tau,1);
  sols(:,i+1)=x(1:2);
end
%
% Plot the solutions.
%
figure(3);
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
print -deps pdafffig3.eps
%eps2pdf('pdafffig3.eps');
