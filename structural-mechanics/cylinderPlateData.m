% Data for the thin cylinder plate problem
% Vinh Phu Nguyen, March 2013
% Cardiff University

addpath ../nurbs-geopdes/inst/
addpath ../post-processing/
addpath ../

Ro   = 2;
Ri   = 1;

%% control points

controlPts = zeros(4,9,2);

% first zKnot 

controlPts(1:3,1,1) = [Ro;0;0];
controlPts(1:3,2,1) = [Ro;Ro;0];
controlPts(1:3,3,1) = [0;Ro;0];
controlPts(1:3,4,1) = [-Ro;Ro;0];
controlPts(1:3,5,1) = [-Ro;0;0];
controlPts(1:3,6,1) = [-Ro;-Ro;0];
controlPts(1:3,7,1) = [0;-Ro;0];
controlPts(1:3,8,1) = [Ro;-Ro;0];
controlPts(1:3,9,1) = [Ro;0;0];

% second zknot

controlPts(1:3,1,2) = [Ri;0;0];
controlPts(1:3,2,2) = [Ri;Ri;0];
controlPts(1:3,3,2) = [0;Ri;0];
controlPts(1:3,4,2) = [-Ri;Ri;0];
controlPts(1:3,5,2) = [-Ri;0;0];
controlPts(1:3,6,2) = [-Ri;-Ri;0];
controlPts(1:3,7,2) = [0;-Ri;0];
controlPts(1:3,8,2) = [Ri;-Ri;0];
controlPts(1:3,9,2) = [Ri;0;0];

controlPts(4,:,:)   = 1;

fac                 = 1/sqrt(2);

controlPts(4,2,1) = fac;
controlPts(4,4,1) = fac;
controlPts(4,6,1) = fac;
controlPts(4,8,1) = fac;
controlPts(4,2,2) = fac;
controlPts(4,4,2) = fac;
controlPts(4,6,2) = fac;
controlPts(4,8,2) = fac;

% homogenous coordinates (x*w,y*w,z*w)

controlPts(1:3,2,1) = controlPts(1:3,2,1)*fac;
controlPts(1:3,4,1) = controlPts(1:3,4,1)*fac;
controlPts(1:3,6,1) = controlPts(1:3,6,1)*fac;
controlPts(1:3,8,1) = controlPts(1:3,8,1)*fac;
controlPts(1:3,2,2) = controlPts(1:3,2,2)*fac;
controlPts(1:3,4,2) = controlPts(1:3,4,2)*fac;
controlPts(1:3,6,2) = controlPts(1:3,6,2)*fac;
controlPts(1:3,8,2) = controlPts(1:3,8,2)*fac;

%% knot vectors 4x1 mesh

uKnot = [0 0 0 1 1 2 2 3 3 4 4 4];
vKnot = [0 0 1 1];

uKnot = uKnot/max(uKnot);

%% build NURBS object

solid = nrbmak(controlPts,{uKnot vKnot});
    
%% k-refinement

% p-refinement to raise order

solid = nrbdegelev(solid,[1 2]); 

% then h-refinement 

refineCount = 2;

for i=1:refineCount
    uKnotVectorU = unique(uKnot);
    uKnotVectorW = unique(vKnot);
    
    % new knots along two directions (uniform)
    
    newKnotsX = uKnotVectorU(1:end-1) + 0.5*diff(uKnotVectorU);
    newKnotsZ = uKnotVectorW(1:end-1) + 0.5*diff(uKnotVectorW);
    newKnots  = {newKnotsX newKnotsZ};
    
    % h-refinement
    
    solid     = nrbkntins(solid,newKnots);    
    uKnot     = cell2mat(solid.knots(1));
    vKnot     = cell2mat(solid.knots(2));
end

% crv = nrbcirc(R,[],deg2rad(0),deg2rad(180));
% nrbplot(crv,80);
% crv = nrbcirc(R-t,[],deg2rad(0),deg2rad(180));
% nrbplot(crv,80);

%%%%%%%%
%% convert NURBS data back to our data structure for analysis

convert2DNurbs

%% Material properties

E  = 3e6;
nu = 0.3;
t  = 0.01; % thickness

%% Boundary condition

q0       = -1;
clamped  = 0;
