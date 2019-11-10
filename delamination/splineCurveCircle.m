addpath ../C_files/
addpath ../fem_util/
addpath ../examples/
addpath ../meshing/

clear all


% initial curve

knotVec    = [0 0 0 1 1 2 2 3 3 4 4 4];
knotVec    = knotVec/max(knotVec);
controlPts = [1 0; 1.4 1.4; 0 1; -1.4 1.4;-1 0;-1.4 -1.4;0 -1; 1.4 -1.4; 1 0];
p          = 2;

weights    = ones(1,9);
weights(2) = 1/sqrt(2);
weights(4) = 1/sqrt(2);
weights(6) = 1/sqrt(2);
weights(8) = 1/sqrt(2);

noCtrPts    = size(controlPts,1);
cp          = zeros(4,noCtrPts);
cp(1:2,:)   = controlPts';
cp(4,:)     = weights;
cp(1,:)   = cp(1,:).*weights;
cp(2,:)   = cp(2,:).*weights;

originalCurve = nrbmak(cp,knotVec);


figure
hold on
nrbctrlplot(originalCurve);
axis equal
axis off

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,'splinecurve.eps',opts)

thickness = 0.15;
% backtracking parameters
alpha = 0.1;
beta  = 0.7;

eps1   = 1e-2;
eps2   = 1e-2;
maxIter = 40;

offsetCurve  = offsetCurveNormal(originalCurve,thickness,alpha,beta,eps2,maxIter);
offsetCurve2 = offsetCurveNormal(originalCurve,2*thickness,alpha,beta,eps2,maxIter);
offsetCurve3 = offsetCurveNormal(originalCurve,3*thickness,alpha,beta,eps2,maxIter);


%% plot result

figure
hold on
%plot(offsetPts(:,1),offsetPts(:,2),'r--','LineWidth',1.9);
%plot(samPts(:,1),   samPts(:,2),'cy--','LineWidth',1.3);

nrbplot(offsetCurve,80);
nrbplot(offsetCurve2,80);
nrbplot(offsetCurve3,80);
nrbplot(originalCurve,80);
nrbctrlplot(originalCurve);
axis equal
axis off


%% 
% extruded surface

srf1 = nrbextrude(originalCurve, [0,1]);
srf2 = nrbextrude(offsetCurve, [0,1]);

figure
hold on
nrbctrlplot(srf1);
nrbctrlplot(srf2);



