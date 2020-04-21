%
% Load the second example problem.  
%
load ex2.mat
x=x0;
y=y0;
z=z0;
taus=logspace(0,-5);
sols=zeros(2,50);
for i=1:50
  [x,y,z]=centpath(A,b,c,x,y,z,taus(i));
  sols(:,i)=x(1:2);
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
xlabel('x_1');
ylabel('x_2');
print -deps centpathfig.eps
%eps2pdf('centpathfig.eps');
