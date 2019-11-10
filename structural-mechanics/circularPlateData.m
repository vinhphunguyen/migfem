%% Geometry data
% circular plate, control points and weights from AV Vuong isogat paper

% plate dimensions
r = 2;

% knots
uKnot = [0 0 0 1 1 1];
vKnot = [0 0 0 1 1 1];

% control points for r=0.5;
controlPts          = zeros(4,3,3);

controlPts(1:2,1,1) = [-sqrt(2)/4 sqrt(2)/4];
controlPts(1:2,2,1) = [-sqrt(2)/2 0];
controlPts(1:2,3,1) = [-sqrt(2)/4 -sqrt(2)/4];

controlPts(1:2,1,2) = [0 sqrt(2)/2];
controlPts(1:2,2,2) = [0;0];
controlPts(1:2,3,2) = [0 -sqrt(2)/2];

controlPts(1:2,1,3) = [sqrt(2)/4 sqrt(2)/4];
controlPts(1:2,2,3) = [sqrt(2)/2 0];
controlPts(1:2,3,3) = [sqrt(2)/4 -sqrt(2)/4];

controlPts(1:2,:,:) = r/0.5*controlPts(1:2,:,:);


% weights
controlPts(4,:,:)   = 1;

fac = sqrt(2)/2;

controlPts(4,2,1) = fac;
controlPts(4,1,2) = fac;
controlPts(4,3,2) = fac;
controlPts(4,2,3) = fac;

controlPts(1:2,2,1) = fac*controlPts(1:2,2,1);
controlPts(1:2,1,2) = fac*controlPts(1:2,1,2);
controlPts(1:2,3,2) = fac*controlPts(1:2,3,2);
controlPts(1:2,2,3) = fac*controlPts(1:2,2,3);

%% build NURBS object

solid = nrbmak(controlPts,{uKnot vKnot});

%% p-refinment

solid = nrbdegelev(solid,[2 2]); % to cubic-cubic NURBS

%% h-refinement

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

nrbkntplot(solid)
nrbctrlplot(solid)
view([0 90])
%% 

convert2DNurbs

res = 200; % resolution for plotting NURBS
plotMesh (controlPts,weights,uKnot,vKnot,p,q,res,'b-','try.eps');

%% Material properties

E  = 30e6;
nu = 0.2;
t  = 0.02; % thickness

%% Boundary condition

q0  = -1.;  % distributed force

clamped = 1;



