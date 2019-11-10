%% Geometry data
%  1/4 model

% plate dimensions
a = 50.0;
b = 50.0; 

% knots
uKnot = [0 0 1 1];
vKnot = [0 0 1 1];

% control points
controlPts          = zeros(4,2,2);

controlPts(1:2,1,1) = [0;0];
controlPts(1:2,2,1) = [a;0;];

controlPts(1:2,1,2) = [0;b];
controlPts(1:2,2,2) = [a;b];

% weights
controlPts(4,:,:)   = 1;

%% build NURBS object

solid = nrbmak(controlPts,{uKnot vKnot});

%% p-refinment

solid = nrbdegelev(solid,[2 2]); % to cubic-cubic NURBS

%% h-refinement

refineLevel = 3;
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

nrbkntplot(solid)
nrbctrlplot(solid)
%% 

convert2DNurbs

res = 200; % resolution for plotting NURBS
plotMesh (controlPts,weights,uKnot,vKnot,p,q,res,'b-','try.eps');

%% Material properties

E  = 1e7;
nu = 0.3;
t  = a/100; % thickness

%% Boundary condition

q0  = -1.;  % distributed force

clamped = 1;



