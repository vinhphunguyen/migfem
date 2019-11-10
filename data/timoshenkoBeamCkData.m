%% Geometry data

% plate dimensions
a = 48.0;
b = 12.0; 

% knots
uKnot = [0 0 1 1];
vKnot = [0 0 1 1];

% control points
controlPts          = zeros(4,2,2);

controlPts(1:2,1,1) = [0;-b/2];
controlPts(1:2,2,1) = [a;-b/2];

controlPts(1:2,1,2) = [0;b/2];
controlPts(1:2,2,2) = [a;b/2];

% weights
controlPts(4,:,:)   = 1;

%% build NURBS object

solid = nrbmak(controlPts,{uKnot vKnot});

%% p-refinment

solid = nrbdegelev(solid,[0 0]); % to cubic-linear NURBS

newKnots  = {[] [0.5]};
solid     = nrbkntins(solid,newKnots);

%% h-refinement

refineLevel = 4;
for i=1:refineLevel
    uKnotVectorU = unique(uKnot);
    uKnotVectorV = unique(vKnot);
    
    % new knots along two directions
    
    newKnotsX = uKnotVectorU(1:end-1) + 0.5*diff(uKnotVectorU);
    newKnotsY = uKnotVectorV(1:end-1) + 0.5*diff(uKnotVectorV);
    
    newKnots  = {newKnotsX []};
    solid     = nrbkntins(solid,newKnots);
    uKnot      = cell2mat(solid.knots(1));
    vKnot      = cell2mat(solid.knots(2));
end

%% 

convert2DNurbs

plotMesh (controlPts,weights,uKnot,vKnot,p,q,50,'b-','try.eps');

