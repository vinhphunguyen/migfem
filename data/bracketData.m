addpath nurbs-geopdes/inst/

%% geometry data
r   = 1;
fac = 1.414213562373095/2;
a   = fac*r;
L   = 8;
b   = L*fac;

%% NURBS data
controlPts          = zeros(4,3,3);
controlPts(1:2,1,1) = [a;a];
controlPts(1:2,2,1) = [0;2*a];
controlPts(1:2,3,1) = [-a;a];

controlPts(1:2,1,2) = [fac*b/2;fac*b/2];
controlPts(1:2,2,2) = [0;fac*b];
controlPts(1:2,3,2) = [-fac*b/2;fac*b/2];

controlPts(1:2,1,3) = [L/2;L/2];
controlPts(1:2,2,3) = [0;L/2];
controlPts(1:2,3,3) = [-L/2;L/2];

controlPts(4,:,:)   = 1;
controlPts(4,2,1)   = fac;
controlPts(1:2,2,1) = controlPts(1:2,2,1)*fac; % homogeneous coord.

% knot vectors
uKnot = [0 0 0 1 1 1];
vKnot = [0 0 0 1 1 1];

%% build NURBS object

% rotation matrix
trans  = vecrotz(pi/2);

trans =[
    0.0000   -1.0000         0         0
    1.0000    0.0000         0         0
    0         0    1.0000         0
    0         0         0    1.0000];

% build first quarter
solid1 = nrbmak(controlPts,{uKnot vKnot});
% then, obtain the remaining by transformation
solid2 = nrbtform(solid1, trans);
solid3 = nrbtform(solid2, trans);
solid4 = nrbtform(solid3, trans);

% plot for visualization
figure
hold on
nrbplot(solid1,[10 10])
nrbplot(solid2,[20 20])
nrbplot(solid3,[20 20])
nrbplot(solid4,[20 20])
view([0 90])
axis equal

%% h-refinement

refineCount = 3;

solid1 = hRefineNURBS(solid1,refineCount);
solid2 = hRefineNURBS(solid2,refineCount);
solid3 = hRefineNURBS(solid3,refineCount);
solid4 = hRefineNURBS(solid4,refineCount);

%% converting to patch data structure

patch1 = convert2DNurbsToPatch(solid1);
patch2 = convert2DNurbsToPatch(solid2);
patch3 = convert2DNurbsToPatch(solid3);
patch4 = convert2DNurbsToPatch(solid4);

noPatches  = 4;
patches(1) = patch1;
patches(2) = patch2;
patches(3) = patch3;
patches(4) = patch4;

%% build mesh structure

nodePattern  = zeros(patch1.noPtsY,patch1.noPtsX);

count = 1;

for i=1:patch1.noPtsY
    for j=1:patch1.noPtsX
        nodePattern(i,j) = count;
        count = count + 1;
    end
end

nodePattern2      = zeros(patch2.noPtsY,patch2.noPtsX);
index             = nodePattern(:,end);
nodePattern2(:,1) = index;
nodeId            = index(end); % maximum number of node of patch1

count = 1;

for i=1:patch2.noPtsY
    for j=2:patch2.noPtsX
        nodePattern2(i,j) = count + nodeId;
        count = count + 1;
    end
end

nodePattern3      = zeros(patch3.noPtsY,patch3.noPtsX);
index             = nodePattern2(:,end);
nodePattern3(:,1) = index;
nodeId            = index(end); % maximum number of node of patch2

count = 1;

for i=1:patch3.noPtsY
    for j=2:patch3.noPtsX
        nodePattern3(i,j) = count + nodeId;
        count = count + 1;
    end
end

nodePattern4        = zeros(patch4.noPtsY,patch4.noPtsX);
index               = nodePattern3(:,end);
index1              = nodePattern (:,1);
nodePattern4(:,1)   = index;
nodePattern4(:,end) = index1;
nodeId              = index(end); % maximum number of node of patch2

count = 1;

for i=1:patch4.noPtsY
    for j=2:patch4.noPtsX-1
        nodePattern4(i,j) = count + nodeId;
        count = count + 1;
    end
end

generateMesh(patch1,nodePattern, nodePattern);
generateMesh(patch2,nodePattern2,nodePattern);
generateMesh(patch3,nodePattern3,nodePattern);
generateMesh(patch4,nodePattern4,nodePattern);

%% plot mesh
figure

for ip=1:4
    uKnot      = patches(ip).uKnot;
    vKnot      = patches(ip).vKnot;
    controlPts = patches(ip).controlPts;
    weights    = patches(ip).weights;
    p          = patches(ip).p;
    q          = patches(ip).q;
    plotMesh (controlPts,weights, uKnot,vKnot,...
        p,q,80, 'b-', 'tem')
end

%% debug only
% vtuFile0 = 'bracket';
%
% % Loop over patches
%
% for ip=1:noPatches
%     index      = patches(ip).index;
%     elRangeU   = patches(ip).elRangeU;
%     elRangeV   = patches(ip).elRangeV;
%     element    = patches(ip).element;
%     elementL   = patches(ip).elementLocal;
%     uKnot      = patches(ip).uKnot;
%     vKnot      = patches(ip).vKnot;
%     controlPts = patches(ip).controlPts;
%     weights    = patches(ip).weights;
%     noPtsX     = patches(ip).noPtsX;
%     noPtsY     = patches(ip).noPtsY;
%     noElems    = patches(ip).noElemsU * patches(ip).noElemsV;
%
%     vtuFile    = strcat(vtuFile0,num2str(ip));
%
%     buildVisualizationMesh;
%
%     sigmaXX = zeros(size(node,1),1);
%     sigmaYY = zeros(size(node,1),1);
%     sigmaXY = zeros(size(node,1),1);
%
%     dispX = zeros(size(node,1),1);
%     dispY = zeros(size(node,1),1);
%
%     VTKPostProcess(node,elementV,2,'Quad4',vtuFile,...
%         [sigmaXX sigmaYY sigmaXY],[dispX dispY]);
% end

