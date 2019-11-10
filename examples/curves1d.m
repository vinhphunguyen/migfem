addpath C_files/
addpath('~/code/xfem-efg-matlab/fem_util');

clear all

knots      = [0 0 0 1 1 1];
controlPts = [0 0; 0 1; 1 1];
p          = 2;

weights1   = [1 1      1]; % b-spline curves
weights2   = [1 0.7071 1]; % NURBS curves

noPts      = 60;
xi         = linspace(0,1,noPts);
sCurve     = zeros(2,noPts);
nCurve     = zeros(2,noPts);


for i=1:noPts
  sCurve(1,i) = NURBSinterpolation(xi(i), p, knots, controlPts(:,1), weights1);
  sCurve(2,i) = NURBSinterpolation(xi(i), p, knots, controlPts(:,2), weights1);
  nCurve(1,i) = NURBSinterpolation(xi(i), p, knots, controlPts(:,1), weights2);
  nCurve(2,i) = NURBSinterpolation(xi(i), p, knots, controlPts(:,2), weights2);
end

hold on
plot(controlPts(:,1),controlPts(:,2),'r-o',...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',10);
plot(sCurve(1,:),sCurve(2,:),'r-','LineWidth',1.8);
plot(nCurve(1,:),nCurve(2,:),'b-','LineWidth',1.8);
axis([0 1 0 1])
axis equal
axis tight

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
exportfig(gcf,'curve1d.eps',opts)



