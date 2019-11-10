a = 0.3; % inner radius
b = 0.5; % outer radius

res = 200; % resolution for plotting NURBS

uKnot = [0 0 0 1 1 1];
vKnot = [0 0 0 1 1 1];

controlPts          = zeros(4,3,3);

controlPts(1:2,1,1) = [a;0];
controlPts(1:2,2,1) = [a;a;];
controlPts(1:2,3,1) = [0;a];

controlPts(1:2,1,2) = [0.5*(a+b);0];
controlPts(1:2,2,2) = [0.5*(a+b); 0.5*(a+b)];
controlPts(1:2,3,2) = [0;0.5*(a+b)];

controlPts(1:2,1,3) = [b;0];
controlPts(1:2,2,3) = [b;b];
controlPts(1:2,3,3) = [0;b];

controlPts(4,:,:)   = 1;

fac                 = 1/sqrt(2);

controlPts(4,2,1) = fac;
controlPts(4,2,2) = fac;
controlPts(4,2,3) = fac;

% homogenous coordinates (x*w,y*w,z*w)

controlPts(1:2,2,1) = controlPts(1:2,2,1)*fac;
controlPts(1:2,2,2) = controlPts(1:2,2,2)*fac;
controlPts(1:2,2,3) = controlPts(1:2,2,3)*fac;

%% build NURBS object

solid = nrbmak(controlPts,{uKnot vKnot});

figure
hold on
nrbplot(solid,[40 40])

%% refinement

refineLevel = 4;
for i=1:refineLevel
    uKnotVectorU = unique(uKnot);
    uKnotVectorV = unique(vKnot);
    
    % new knots along two directions
    
    newKnotsX = uKnotVectorU(1:end-1) + 0.5*diff(uKnotVectorU);
    newKnotsY = uKnotVectorV(1:end-1) + 0.5*diff(uKnotVectorV);
    
    newKnots  = {newKnotsX newKnotsY};
    solid     = nrbkntins(solid,newKnots);
    uKnot      = cell2mat(solid.knots(1));
    vKnot      = cell2mat(solid.knots(2));
end

%%

% newKnotsX = 0.5;
%     newKnotsY = 0.5;
%
%     newKnots  = {newKnotsX newKnotsY};
%     solid     = nrbkntins(solid,newKnots);

convert2DNurbs

plotMesh (controlPts,weights,uKnot,vKnot,p,q,res,'r-','try.eps');

%% Bezier extraction operators

[C,Cxi,Cet]  = bezierExtraction2D(uKnot,vKnot,p,q);

noCtrPts       = noPtsX * noPtsY;
noDofs         = noCtrPts * 2;


%%

% generateIGA2DMesh
% bezierPoints = [];
% 
% for e=1:noElems
%     sctr   = element(e,:);         %  element scatter vector
%     
%     Ce     = C(:,:,e);             % element Bezier extraction operator
%     we     = diag(weights(sctr,:));% element weights
%     pts    = controlPts(sctr,:);   % element nodes
%     Wb     = Ce'*weights(sctr,:);  % element Bezier weights
%     
%     bezierPts = inv(diag(Wb))*Ce'*we*pts;
%     bezierPoints = [bezierPoints;bezierPts];
% end
% 
% % remove sharing control points
% 
% bezierPoints = unique_no_sort_rows(bezierPoints);
% 
% bezierElem=[1 4 7 2 5 8 3 6 9;
%             7 16 19 8 17 20 9 18 21;
%             3 6 9 10 12 14 11 13 15;
%             9 18 21 14 22 24 15 23 25];
% 
% plot_mesh(bezierPoints,bezierElem,'Q9','b-',1);
% 
% plot(bezierPoints(:,1),bezierPoints(:,2),'rs','MarkerSize',13,'MarkerFaceColor','r');

%%
