% An example of the de Casteljau algoritm

addpath C_files/
addpath('~/code/xfem-efg-matlab/fem_util');

clear all

knots      = [0 0 0 1 1 1];
controlPts = [1 1; 
              3 5; 
              5 2];
p          = 2;

weights  = [1 1 1]; % b-spline curves


noPts      = 60;
xi         = linspace(0,1,noPts);
sCurve     = zeros(2,noPts);


for i=1:noPts
  sCurve(1,i) = NURBSinterpolation(xi(i), p, knots, controlPts(:,1), weights);
  sCurve(2,i) = NURBSinterpolation(xi(i), p, knots, controlPts(:,2), weights);
end

hold on
plot(controlPts(:,1),controlPts(:,2),'b-o',...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',10, 'LineWidth',1.2);
plot(sCurve(1,:),sCurve(2,:),'r-','LineWidth',1.8);

%axis([0 1 0 1])
axis equal
axis tight
axis off


% Compute point at t=0.5 using the 
% de Casteljau algorithm

t = 0.4;

p0 = controlPts(1,:);
p1 = controlPts(2,:);
p2 = controlPts(3,:);

p01 = (1-t)*p0 + (t)*p1;
p11 = (1-t)*p1 + (t)*p2;

p02 = (1-t)*p01 + (t)*p11;

p0111 = [p01;p11];

plot(p0111(:,1),p0111(:,2),'r-s',...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',10);

            
plot(p02(:,1),p02(:,2),'bs',...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','r',...
                'MarkerSize',10);

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
exportfig(gcf,'curve1d.eps',opts)



