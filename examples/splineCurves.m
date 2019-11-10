addpath ../C_files/
addpath ../fem_util/

clear all

knots      = [0 0 0 1 2 3 4 5 5 5];
knots1     = [0 0 0 1 2 3 4 4 5 5 5];

uniqueKnot = unique(knots1);

controlPts  = [0 1; 0.5 3; 1.5 3; 2.5 0; 3.5 2; 4.5 1; 5 4];
controlPts1 = [0 1; 0.5 3; 1.5 3; 2.5 0; 3.5 2; 4. 0.5; 4.5 1;5 4];

p          = 2;

weights1   = [1 1 1 1 1 1 1]; % b-spline curves
weights2   = [1 1 1 1 1 1 1 1]; % b-spline curves


noPts      = 80;
xi         = linspace(0,max(knots),noPts);
xi1        = linspace(0,max(knots1),noPts);
sCurve     = zeros(2,noPts);
nCurve     = zeros(2,noPts);
iKnot      = zeros(2,length(uniqueKnot));


for i=1:noPts
  sCurve(1,i) = NURBSinterpolation(xi(i), p, knots, controlPts(:,1), weights1);
  sCurve(2,i) = NURBSinterpolation(xi(i), p, knots, controlPts(:,2), weights1);
  nCurve(1,i) = NURBSinterpolation(xi1(i), p, knots1, controlPts1(:,1), weights2);
  nCurve(2,i) = NURBSinterpolation(xi1(i), p, knots1, controlPts1(:,2), weights2);
end

for i=1:length(uniqueKnot)
  iKnot(1,i) = NURBSinterpolation(uniqueKnot(i), p, knots1, controlPts1(:,1), weights2);
  iKnot(2,i) = NURBSinterpolation(uniqueKnot(i), p, knots1, controlPts1(:,2), weights2);
end

hold on
plot(controlPts(:,1),controlPts(:,2),'r-o',...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',10);
plot(sCurve(1,:),sCurve(2,:),'b-','LineWidth',1.8);
%plot(nCurve(1,:),nCurve(2,:),'b-','LineWidth',1.8);
axis([0 1 0 1])
axis equal
axis tight

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
exportfig(gcf,'splinecurve.eps',opts)

figure
hold on
plot(controlPts1(:,1),controlPts1(:,2),'r-o',...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',10);
plot(nCurve(1,:),nCurve(2,:),'b-','LineWidth',1.8);
plot(iKnot(1,:),iKnot(2,:),'ro','LineWidth',1.8);
axis([0 1 0 1])
axis equal
axis tight

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
exportfig(gcf,'splinecurve1.eps',opts)


