% shell-solid coupling problem
% solid = 1 trivariate NURBS patch
% shell = 1 bivariate NURBS patch
% there are common control points at their intersection

addpath ../fem_util/
addpath ../nurbs-geopdes/inst/
addpath ../nurbs-util/
addpath ../meshing/
addpath ../fem-functions/
addpath ../post-processing/
addpath ../xiga/
addpath ../nurbs-geopdes/inst/
addpath ../delamination/
addpath ../structural-mechanics/
addpath ../xml_toolbox/
addpath ../C_files/

clear all

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solid patch

a = 10;
b = 2;
c = 10;

controlPts=zeros(4,2,2,2);
controlPts(1:3,1,1,1)=[0;0;0];
controlPts(1:3,2,1,1)=[a;0;0];
controlPts(1:3,1,2,1)=[0;b;0];
controlPts(1:3,2,2,1)=[a;b;0];
controlPts(1:3,1,1,2)=[0;0;c];
controlPts(1:3,2,1,2)=[a;0;c];
controlPts(1:3,1,2,2)=[0;b;c];
controlPts(1:3,2,2,2)=[a;b;c];
controlPts(4,:,:)   = 1;
uKnot = [0 0 1 1];
vKnot = [0 0 1 1];
wKnot = [0 0 1 1];

solid = nrbmak(controlPts,{uKnot vKnot wKnot});

uKnot     = cell2mat(solid.knots(1));
vKnot     = cell2mat(solid.knots(2));
wKnot     = cell2mat(solid.knots(3));

refineCountX = 1;
for i=1:refineCountX
    uKnotVectorU = unique(uKnot);
    uKnotVectorV = unique(vKnot);
    uKnotVectorW = unique(wKnot);
    
    % new knots along two directions (uniform)
    
    newKnotsX = uKnotVectorU(1:end-1) + 0.5*diff(uKnotVectorU);
    newKnotsY = uKnotVectorV(1:end-1) + 0.5*diff(uKnotVectorV);
    newKnotsZ = uKnotVectorW(1:end-1) + 0.5*diff(uKnotVectorW);
    newKnots  = {newKnotsX [] newKnotsZ};
    
    % h-refinement
    solid     = nrbkntins(solid,newKnots);
    uKnot     = cell2mat(solid.knots(1));
    vKnot     = cell2mat(solid.knots(2));
    wKnot     = cell2mat(solid.knots(3));
end
solidMesh  = buildIGA3DMesh(solid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% shell patch
l = 5;
controlPts=zeros(4,2,2);
controlPts(1:3,1,1)=[0.5*a;b;0];
controlPts(1:3,2,1)=[0.5*a;b;c];
controlPts(1:3,1,2)=[0.5*a;b+l;0];
controlPts(1:3,2,2)=[0.5*a;b+l;c];
controlPts(4,:,:)   = 1;
uKnot = [0 0 1 1];
vKnot = [0 0 1 1];

shell = nrbmak(controlPts,{uKnot vKnot});

figure
hold on
nrbplot(solid,[10 1 10]);
nrbctrlplot(shell);
axis equal
axis off



