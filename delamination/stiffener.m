addpath ../C_files/
addpath ../fem_util/
addpath ../examples/
addpath ../meshing/
addpath ../nurbs-util/

clear all
eps=1e-8;
refCount = [3 2 2];


% initial curve (quadratic)

knotVec     = [0 0 0 1 1 1];
controlPts  = [0 0; 1 -1 ;2 0];
weights     = [1 1 1]; % b-spline curves
p           = 2;

noCtrPts    = size(controlPts,1);
cp          = zeros(4,noCtrPts);
cp(1:2,:)   = controlPts';
cp(4,:)     = weights;

skinCurve = nrbmak(cp,knotVec);


figure
hold on
nrbctrlplot(skinCurve);
axis equal
axis off

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,'splinecurve.eps',opts)

%% insert knots to have compatibility

knots = [0.3 0.3 0.4 0.4 0.7 0.7 0.6 0.6];

skinCurve     = nrbkntins(skinCurve,knots);

figure
hold on
nrbctrlplot(skinCurve);
axis equal
axis off

p4 = skinCurve.coefs(1:2,4)';
p3 = skinCurve.coefs(1:2,3)';
p5 = skinCurve.coefs(1:2,5)';
p7 = skinCurve.coefs(1:2,7)';
p8 = skinCurve.coefs(1:2,8)';
p9 = skinCurve.coefs(1:2,9)';

%% stiffener curve

stiffKnots  = [0 0 0 1 1 3 3 4 4 4];
stiffCtrPts = [p3; p4; p5; 1 0.2; p7; p8; p9];
stiffWgts   = [1 1 1 1 1 1 1]; % b-spline curves
stiffOrder  = 2;

noCtrPts    = size(stiffCtrPts,1);
cp          = zeros(4,noCtrPts);
cp(1:2,:)   = stiffCtrPts';
cp(4,:)     = stiffWgts;

stiffCurve = nrbmak(cp,stiffKnots);


figure
hold on
nrbctrlplot(stiffCurve);
nrbctrlplot(skinCurve);
axis equal
axis off

%% offsetting

% backtracking parameters
alpha = 0.1;
beta  = 0.7;

eps1   = 1e-2;
eps2   = 1e-2;

maxIter = 100;

skinThickness  = -0.1;
stiffThickness = 0.05;

[offsetSkinCurve,offsetPts]= offsetCurveNormal(skinCurve,skinThickness,alpha,beta,eps1,maxIter);
offsetStiffCurve = offsetCurve(stiffCurve,stiffThickness,alpha,beta,eps2,maxIter);


figure
hold on
nrbctrlplot(skinCurve);
%nrbctrlplot(stiffCurve);
nrbctrlplot(offsetSkinCurve);
nrbctrlplot(offsetStiffCurve);
axis equal
axis off

pause

%% make surfaces from curves

skin = surfaceFromTwoCurves(skinCurve, offsetSkinCurve);
stif = surfaceFromTwoCurves(stiffCurve,offsetStiffCurve);

figure
hold on
nrbplot(skin,[20 2])
nrbplot(stif,[20 2])
axis off
view([0 90])

%% make volume from surfaces

skinVol = nrbextrude(skin, [0,0,1]);
stifVol = nrbextrude(stif, [0,0,1]);

figure
hold on
nrbplot(skinVol,[20 2 2])
nrbplot(stifVol,[20 2 2])
axis off

% k-refinement
ods1 = skinVol.order -1;
ods2 = stifVol.order -1;
skinVol = doKRefinementSolid(skinVol,ods1,refCount);
stifVol = doKRefinementSolid(stifVol,ods2,refCount);

% FE mesh from NURBS solids
skinMesh3D = buildIGA3DMesh(skinVol);
stifMesh3D = buildIGA3DMesh(stifVol);
% visualization meshes (one for 1 patch)
vMeshSkin3D = buildVisualizationMesh3D(skinVol);
vMeshStif3D = buildVisualizationMesh3D(stifVol);

figure; hold on;
plot_mesh(vMeshSkin3D.node,vMeshSkin3D.element,'B8','b-',1.1);
plot_mesh(vMeshStif3D.node,vMeshStif3D.element,'B8','r-',1.1);

stifMesh3D.globElems = stifMesh3D.globElems + size(skinMesh3D.controlPts,1);

% store the meshes in data structure 'mesh'
mesh{1} = skinMesh3D;
mesh{2} = stifMesh3D;

vmesh{1} = vMeshSkin3D;
vmesh{2} = vMeshStif3D;

%Dirichlet nodes
% must in column vector

pp1 = find(abs(skinMesh3D.controlPts(:,1)-offsetPts(1,1))  <eps);
pp2 = find(abs(skinMesh3D.controlPts(:,1)-offsetPts(end,1))<eps);

xnodes = [pp1;pp2];
ynodes = xnodes;
znodes = xnodes;

aa1 = find(abs(skinMesh3D.controlPts(:,1)-p3(1))<eps);
aa2 = find(abs(skinMesh3D.controlPts(:,1)-p5(1))<eps);
aa3 = find(abs(skinMesh3D.controlPts(:,1)-p7(1))<eps);
aa4 = find(abs(skinMesh3D.controlPts(:,1)-p9(1))<eps);

bb1 = find(abs(stifMesh3D.controlPts(:,1)-p3(1))<eps);
bb2 = find(abs(stifMesh3D.controlPts(:,1)-p5(1))<eps);
bb3 = find(abs(stifMesh3D.controlPts(:,1)-p7(1))<eps);
bb4 = find(abs(stifMesh3D.controlPts(:,1)-p9(1))<eps);

aa  = [];
bb  = [];

for i=1:length(aa1)
    s = aa1(i);
    e = aa2(i);
    aa = [aa; s:e];
    s = bb1(i);
    e = bb2(i);
    bb = [bb; s:e];
end

for i=1:length(aa3)
    s = aa3(i);
    e = aa4(i);
    aa = [aa; s:e];
    s = bb3(i);
    e = bb4(i);
    bb = [bb; s:e];
end

bb = bb + size(skinMesh3D.controlPts,1);

% datastructure, 'data', contains the whole mesh

data.mesh     = mesh;
data.vmesh    = vmesh;
data.pntCount = size(mesh{1}.controlPts,1) + size(mesh{2}.controlPts,1);
data.xnodes   = xnodes;
data.ynodes   = ynodes;
data.znodes   = znodes;


%%

% k refinement here

% refCount = 0;
% skin = doKRefinementSurface(skin,skin.order(1)-1,skin.order(2)-1,refCount);
% stif = doKRefinementSurface(stif,stif.order(1)-1,stif.order(2)-1,refCount);

% [skinmesh] = buildIGA2DMeshForSurface (skin);
% [stifmesh] = buildIGA2DMeshForSurface (stif);
%
% vMeshSkin = buildVisualizationMesh2D(skin);
% vMeshStif = buildVisualizationMesh2D(stif);
%
% eps=1e-8;
% aa1 = find(abs(skinmesh.controlPts(:,1)-p3(1))<eps);
% aa2 = find(abs(skinmesh.controlPts(:,1)-p5(1))<eps);
% aa3 = find(abs(skinmesh.controlPts(:,1)-p7(1))<eps);
% aa4 = find(abs(skinmesh.controlPts(:,1)-p9(1))<eps);
%
% bb1 = find(abs(stifmesh.controlPts(:,1)-p3(1))<eps);
% bb2 = find(abs(stifmesh.controlPts(:,1)-p5(1))<eps);
% bb3 = find(abs(stifmesh.controlPts(:,1)-p7(1))<eps);
% bb4 = find(abs(stifmesh.controlPts(:,1)-p9(1))<eps);
%
% aa  = [aa1:aa2 aa3:aa4];
% bb  = [bb1:bb2 bb3:bb4];
%
% figure
% hold on
% n1 = plot(skinmesh.controlPts(aa,1),skinmesh.controlPts(aa,2),'ro' );
% n2 = plot(stifmesh.controlPts(bb,1),stifmesh.controlPts(bb,2),'bs' );
% set(n1,'MarkerSize',12, 'MarkerFaceColor','g','MarkerEdgeColor','k');
% set(n2,'MarkerSize',12, 'MarkerFaceColor','cy','MarkerEdgeColor','k');
% plot_mesh(vMeshSkin.node,vMeshSkin.element,'Q4','b-',1.1);
% plot_mesh(vMeshStif.node,vMeshStif.element,'Q4','r-',1.1);
%
% stifmesh.globElems = stifmesh.globElems + size(skinmesh.controlPts,1);
% bb = bb + + size(skinmesh.controlPts,1);
%
% % for e=1:size(stifmesh.globElems,1)
% %     sctr  = stifmesh.locElems(e,:);
% %     sctrg = stifmesh.globElems(e,:);
% %     [C,IA,IB] = intersect(sctr,bb);
% %     remainIds = setdiff(1:length(sctr),IA);
% %     sctr(IA)  = aa(IB);
% %     sctr(remainIds)  = sctrg(remainIds);
% %     stifmesh.globElems(e,:) = sctr;
% % end
%
%
%
% mesh{1} = skinmesh;
% mesh{2} = stifmesh;
%
% %Dirichlet nodes
% % must in column vector
% xnodes = [size(skinmesh.controlPts,1)-skinmesh.noPtsX+1 size(skinmesh.controlPts,1)]';
% ynodes = xnodes;
%
% % datastructure contains the whole mesh
%
% data.mesh     = mesh;
% data.pntCount = size(mesh{1}.controlPts,1) + size(mesh{2}.controlPts,1);
% data.xnodes   = xnodes;
% data.ynodes   = ynodes;






