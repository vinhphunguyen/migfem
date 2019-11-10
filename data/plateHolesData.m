% data for plate with a hole Ck elements
% h and k-refinement can be used.
% Used for a Finite Cell Method.
%
% Vinh Phu Nguyen
% Johns Hopkins University

%%

addpath('~/code/xfem-efg-matlab/fem_util');
addpath ../nurbs-geopdes/inst/
addpath ../nurbs-util/
addpath ../meshing/
addpath ../fem-functions/
addpath ../post-processing/
addpath ../xiga/

%clear all

L = 4; % half plate width

controlPts          = zeros(4,2,2);

controlPts(1:2,1,1) = [-L;-L];
controlPts(1:2,2,1) = [L;-L];
controlPts(1:2,1,2) = [-L;L];
controlPts(1:2,2,2) = [L;L];

controlPts(4,:,:)   = 1;

uKnot = [0 0 1 1];
vKnot = [0 0 1 1];

%% build NURBS object

solid = nrbmak(controlPts,{uKnot vKnot});

% and evaluate order

solid = nrbdegelev(solid,[2 2]);

% h-refinement

refineCount = 5;

for i=1:refineCount
    uKnotVectorU = unique(uKnot);
    uKnotVectorV = unique(vKnot);
    
    % new knots along two directions (uniform)
    
    newKnotsX = uKnotVectorU(1:end-1) + 0.5*diff(uKnotVectorU);
    newKnotsY = uKnotVectorV(1:end-1) + 0.5*diff(uKnotVectorV);
    newKnots  = {newKnotsX newKnotsY};
    
    % h-refinement
    
    solid     = nrbkntins(solid,newKnots);
    
    uKnot      = cell2mat(solid.knots(1));
    vKnot      = cell2mat(solid.knots(2));
end

%%

convert2DNurbs

noCtrPts       = noPtsX * noPtsY;
noDofs         = noCtrPts * 2;

%% generate element connectivity ...

generateIGA2DMesh

% plot the mesh

buildVisualizationMesh;

% circular hole data

r  = 1.2; % radius
xc = -1.5; % x coord of center
yc = 0; % y coord of center

VOID = [-2 -2 r;
        2 -2 r;
        -2 2 r;
        2 2 r];
    
    %VOID = [0 0 r];

levelSetVoids

% Plot mesh and enriched nodes to check
figure
hold on
plot_mesh(node,elementV,'Q4','b-');
% plot the circle
theta = 0:0.01:2*pi;
xo = xc + r*cos(theta) ;
yo = yc + r*sin(theta) ;
plot(xo,yo,'k-','Linewidth',1.9);
% plot elements cut by the circle
plot_mesh(node,elementV(splitElems,:),'Q4','r-');
plot_mesh(node,elementV(inactiveElems,:),'Q4','c-*');
%plot(gps(:,1),gps(:,2),'+');
