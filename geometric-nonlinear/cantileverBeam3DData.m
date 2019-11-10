%% Geometry data

% plate dimensions
a = 10.0; % length
b = 1.0;  % width 
t = 0.1;  % thickness

% knots
uKnot = [0 0 1 1];
vKnot = [0 0 1 1];
wKnot = [0 0 1 1];

% control points
controlPts          = zeros(4,2,2,2);

controlPts(1:3,1,1,1) = [0;0;0];
controlPts(1:3,2,1,1) = [a;0;0];

controlPts(1:3,1,2,1) = [0;b;0];
controlPts(1:3,2,2,1) = [a;b;0];

controlPts(1:3,1,1,2) = [0;0;t];
controlPts(1:3,2,1,2) = [a;0;t];

controlPts(1:3,1,2,2) = [0;b;t];
controlPts(1:3,2,2,2) = [a;b;t];

% weights
controlPts(4,:,:)   = 1;

%% build NURBS object

solid = nrbmak(controlPts,{uKnot vKnot wKnot});

%% p-refinment

solid = nrbdegelev(solid,[2 0 0]); % to cubic-linear NURBS

%% h-refinement

refineLevel = 4;
for i=1:refineLevel
    uKnotVectorU = unique(uKnot);
    uKnotVectorV = unique(vKnot);
    
    % new knots along two directions
    
    newKnotsX = uKnotVectorU(1:end-1) + 0.5*diff(uKnotVectorU);
    newKnotsY = uKnotVectorV(1:end-1) + 0.5*diff(uKnotVectorV);
    
    newKnots  = {newKnotsX {} {}};
    solid     = nrbkntins(solid,newKnots);
    uKnot      = cell2mat(solid.knots(1));
    vKnot      = cell2mat(solid.knots(2));
end

nrbkntplot(solid)
%% 

convert3DNurbs


