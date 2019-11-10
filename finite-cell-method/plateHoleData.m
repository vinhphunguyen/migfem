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

%clear all

L = 4; % half plate width

controlPts          = zeros(4,2,2);

controlPts(1:2,1,1) = [0;0];
controlPts(1:2,2,1) = [L;0];
controlPts(1:2,1,2) = [0;L];
controlPts(1:2,2,2) = [L;L];

controlPts(4,:,:)   = 1;

uKnot = [0 0 1 1];
vKnot = [0 0 1 1];

%% build NURBS object

solid = nrbmak(controlPts,{uKnot vKnot});

% and evaluate order

solid = nrbdegelev(solid,[5 5]);

% h-refinement

refineCount = 4;

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

r  = 1.7; % radius
xc = 0; % x coord of center
yc = 0; % y coord of center

for i = 1 : length(node)
    x      = node(i,1);
    y      = node(i,2);
    d      = sqrt((x-xc)^2+(y-yc)^2);
    ls(i)  = d - r; % level set
end

count1 = 0;
count2 = 0;

% loop over elements

for iel = 1 : length(elementV)
    sctr    = elementV(iel,:);
    phi     = ls(sctr);
    
    if    ( max(phi)*min(phi) < 0 )
        count1              = count1 + 1 ; % one split element
        splitElems(count1)  = iel;
    elseif max(phi) < 0
        count2                 = count2 + 1 ; % one inactive element
        inactiveElems(count2)  = iel;
    end
end



