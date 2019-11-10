addpath ../C_files/
addpath ../fem_util/
addpath ../examples/

clear all

knotVec     = [0 0 0 1 1 1];
controlPts  = [0 0; 1 4; 2 0;];
weights     = [1 1 1]; % b-spline curves



noCtrPts    = size(controlPts,1);
cp          = zeros(4,noCtrPts);
cp(1:2,:)   = controlPts';
cp(4,:)     = weights;
cp(1,:)   = cp(1,:).*weights;
cp(2,:)   = cp(2,:).*weights;

originalCurve = nrbmak(cp,knotVec);

% backtracking parameters
alpha = 0.1;
beta  = 0.7;

eps1   = 1e-2;
eps2   = 1e-3;
maxIter = 100;

thickness = 0.2;

offsetCurve  = offsetCurveNormal(originalCurve,thickness,alpha,beta,eps2,maxIter);



%% plot result

figure
hold on
plot(offsetPts(:,1),offsetPts(:,2),'r--','LineWidth',1.9);
%plot(samPts(:,1),   samPts(:,2),'cy--','LineWidth',1.3);

nrbplot(offsetCurve,80);
nrbplot(originalCurve,80);
%nrbctrlplot(originalCurve);
axis equal
axis off

figure
hold on
nrbctrlplot(originalCurve);
axis equal
nrbctrlplot(offsetCurve);


