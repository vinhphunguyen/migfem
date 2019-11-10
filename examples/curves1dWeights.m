addpath C_files/
addpath('~/code/xfem-efg-matlab/fem_util');

clear all

knots      = [0 0 0 1 2 3 4 4 4];
controlPts = [0 0; 1 0; 1 1;2 1;2 0; 3 0];
p          = 2;

weights1   = [1 1 1 1 1 1]; % b-spline curves
weights2   = [1 1 1 0.6 1 1]; % NURBS curves
weights3   = [1 1 1 0.3 1 1]; % NURBS curves

noPts      = 60;
xi         = linspace(0,1,noPts);
sCurve1    = zeros(2,noPts);
sCurve2    = zeros(2,noPts);
sCurve3    = zeros(2,noPts);

knots = knots/max(knots);

for i=1:noPts
  sCurve1(1,i) = NURBSinterpolation(xi(i), p, knots, controlPts(:,1), weights1);
  sCurve1(2,i) = NURBSinterpolation(xi(i), p, knots, controlPts(:,2), weights1);
  sCurve2(1,i) = NURBSinterpolation(xi(i), p, knots, controlPts(:,1), weights2);
  sCurve2(2,i) = NURBSinterpolation(xi(i), p, knots, controlPts(:,2), weights2);
  sCurve3(1,i) = NURBSinterpolation(xi(i), p, knots, controlPts(:,1), weights3);
  sCurve3(2,i) = NURBSinterpolation(xi(i), p, knots, controlPts(:,2), weights3);
end

hold on
plot(controlPts(:,1),controlPts(:,2),'r-o',...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',10);
plot(sCurve1(1,:),sCurve1(2,:),'r-','LineWidth',1.8);
plot(sCurve2(1,:),sCurve2(2,:),'b-','LineWidth',1.8);
plot(sCurve3(1,:),sCurve3(2,:),'c-','LineWidth',1.8);
%axis([0 1 0 1])
axis equal
axis tight

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
exportfig(gcf,'curve1d.eps',opts)



