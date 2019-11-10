addpath ../C_files/
addpath ../fem_util/
addpath ../examples/
addpath ../meshing/
addpath ../nurbs-util/

clear all
eps=1e-8;
no = 0; % 3=> 8 plies
refCount = [0 no 0];

height = 4;   % height of stiffener
samLen = 50;  % length of sample
width  = 50;   % width of sample
curv   = -0; % large => more curvature

% initial curve (quadratic)

knotVec     = [0 0 0 1 1 1];
controlPts  = [0 0; samLen/2 curv ;samLen 0];
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

knots = [0.1 0.1 0.2 0.2 0.3 0.3 0.4 0.4 0.6 0.6 0.7 0.7 0.8 0.8 0.9 0.9];

skinCurve     = nrbkntins(skinCurve,knots);

figure
hold on
nrbctrlplot(skinCurve);
axis equal
axis off

p6  = skinCurve.coefs(1:2,6)';
p14 = skinCurve.coefs(1:2,14)';

p3 = skinCurve.coefs(1:2,3)';
p4 = skinCurve.coefs(1:2,4)';
p5 = skinCurve.coefs(1:2,5)';
p7 = skinCurve.coefs(1:2,7)';
p8 = skinCurve.coefs(1:2,8)';
p9 = skinCurve.coefs(1:2,9)';
p11 = skinCurve.coefs(1:2,11)';
p12 = skinCurve.coefs(1:2,12)';
p13 = skinCurve.coefs(1:2,13)';
p15 = skinCurve.coefs(1:2,15)';
p16 = skinCurve.coefs(1:2,16)';
p17 = skinCurve.coefs(1:2,17)';

%% stiffener curves

xii=0.25;
[N dNdxi] = NURBS1DBasisDers(xii,p,knotVec,weights);
s   = findspan(noCtrPts-1,p,xii,knotVec);
pts = controlPts(s-1:s-1+p,:);
x   = N    * pts;
dx  = dNdxi* pts;
dx  = dx/norm(dx);
pnt1 = x + [-dx(2) dx(1)]*height;

xii=0.75;
[N dNdxi] = NURBS1DBasisDers(xii,p,knotVec,weights);
s   = findspan(noCtrPts-1,p,xii,knotVec);
pts = controlPts(s-1:s-1+p,:);
x   = N    * pts;
dx  = dNdxi* pts;
dx  = dx/norm(dx);
pnt2 = x + [-dx(2) dx(1)]*height;

stiffKnots  = [0 0 0 1 1 3 3 4 4 4];
stiffCtrPts = [p3; p4; p5; pnt1; p7; p8; p9];
stiffWgts   = [1 1 1 1 1 1 1]; % b-spline curves
stiffOrder  = 2;

noCtrPts    = size(stiffCtrPts,1);
cp          = zeros(4,noCtrPts);
cp(1:2,:)   = stiffCtrPts';
cp(4,:)     = stiffWgts;

stiffCurve1 = nrbmak(cp,stiffKnots);

stiffKnots  = [0 0 0 1 1 3 3 4 4 4];
stiffCtrPts = [p11; p12; p13; pnt2; p15; p16; p17];
stiffWgts   = [1 1 1 1 1 1 1]; % b-spline curves
stiffOrder  = 2;

noCtrPts    = size(stiffCtrPts,1);
cp          = zeros(4,noCtrPts);
cp(1:2,:)   = stiffCtrPts';
cp(4,:)     = stiffWgts;

stiffCurve2 = nrbmak(cp,stiffKnots);

figure
hold on
nrbctrlplot(stiffCurve1);
nrbctrlplot(stiffCurve2);
nrbctrlplot(skinCurve);
axis equal
axis off

%% offsetting

% backtracking parameters
alpha = 0.1;
beta  = 0.7;

eps1   = 1e-2;
eps2   = 1e-2;

maxIter = 30;

skinThickness  = -0.8;
stiffThickness = 0.8;

[offsetSkinCurve,offsetPts]= offsetCurveNormal(skinCurve,skinThickness,alpha,beta,eps1,maxIter);
offsetStiffCurve1 = offsetCurve(stiffCurve1,stiffThickness,alpha,beta,eps2,maxIter);
offsetStiffCurve2 = offsetCurve(stiffCurve2,stiffThickness,alpha,beta,eps2,maxIter);


figure
hold on
nrbctrlplot(skinCurve);
nrbctrlplot(stiffCurve1);
nrbctrlplot(stiffCurve2);
nrbctrlplot(offsetSkinCurve);
nrbctrlplot(offsetStiffCurve1);
nrbctrlplot(offsetStiffCurve2);
axis equal
axis off

pause

%% make surfaces from curves

skin  = surfaceFromTwoCurves(skinCurve, offsetSkinCurve);
stif1 = surfaceFromTwoCurves(stiffCurve1,offsetStiffCurve1);
stif2 = surfaceFromTwoCurves(stiffCurve2,offsetStiffCurve2);

figure
hold on
nrbplot(skin,[20 2])
nrbplot(stif1,[20 2])
nrbplot(stif2,[20 2])
axis off
view([0 90])

%% make volume from surfaces

skinVol  = nrbextrude(skin,  [0,0,width]);
stifVol1 = nrbextrude(stif1, [0,0,width]);
stifVol2 = nrbextrude(stif2, [0,0,width]);

figure
hold on
nrbplot(skinVol,[20 2 10])
nrbplot(stifVol1,[20 2 10])
nrbplot(stifVol2,[20 2 10])
axis off

% k-refinement
ods1 = [2 1 1];
ods2 = [2 1 1];
ods3 = [2 1 1];

skinVol  = doKRefinementSolid(skinVol, ods1,refCount);
stifVol1 = doKRefinementSolid(stifVol1,ods2,refCount);
stifVol2 = doKRefinementSolid(stifVol2,ods3,refCount);

% FE mesh from NURBS solids
skinMesh3D  = buildIGA3DMesh(skinVol);
stifMesh3D1 = buildIGA3DMesh(stifVol1);
stifMesh3D2 = buildIGA3DMesh(stifVol2);
% visualization meshes (one for 1 patch)
vMeshSkin3D  = buildVisualizationMesh3D(skinVol);
vMeshStif3D1 = buildVisualizationMesh3D(stifVol1);
vMeshStif3D2 = buildVisualizationMesh3D(stifVol2);

figure; hold on;
plot_mesh(vMeshSkin3D.node,vMeshSkin3D.element,'B8','b-',1.1);
plot_mesh(vMeshStif3D1.node,vMeshStif3D1.element,'B8','r-',1.1);
plot_mesh(vMeshStif3D2.node,vMeshStif3D2.element,'B8','r-',1.1);

stifMesh3D1.globElems = stifMesh3D1.globElems + size(skinMesh3D.controlPts,1);
stifMesh3D2.globElems = stifMesh3D2.globElems + size(skinMesh3D.controlPts,1) ...
    + size(stifMesh3D1.controlPts,1);




%Dirichlet nodes
% must in column vector

pp1 = find(abs(skinMesh3D.controlPts(:,1)-offsetPts(1,1))  <eps);
pp2 = find(abs(skinMesh3D.controlPts(:,1)-offsetPts(end,1))<eps);

xnodes = [pp1;pp2];
ynodes = xnodes;
znodes = xnodes;

aa1 = intersect(find(abs(skinMesh3D.controlPts(:,1)-p3(1))<eps),...
                find(abs(skinMesh3D.controlPts(:,2)-p3(2))<eps));
aa2 = intersect(find(abs(skinMesh3D.controlPts(:,1)-p5(1))<eps),...
                find(abs(skinMesh3D.controlPts(:,2)-p5(2))<eps));
aa3 = intersect(find(abs(skinMesh3D.controlPts(:,1)-p7(1))<eps),...
                find(abs(skinMesh3D.controlPts(:,2)-p7(2))<eps));
aa4 = intersect(find(abs(skinMesh3D.controlPts(:,1)-p9(1))<eps),...
                find(abs(skinMesh3D.controlPts(:,2)-p9(2))<eps));
aa5 = intersect(find(abs(skinMesh3D.controlPts(:,1)-p11(1))<eps),...
                find(abs(skinMesh3D.controlPts(:,2)-p11(2))<eps));
aa6 = intersect(find(abs(skinMesh3D.controlPts(:,1)-p13(1))<eps),...
                find(abs(skinMesh3D.controlPts(:,2)-p13(2))<eps));
aa7 = intersect(find(abs(skinMesh3D.controlPts(:,1)-p15(1))<eps),...
                find(abs(skinMesh3D.controlPts(:,2)-p15(2))<eps));
aa8 = intersect(find(abs(skinMesh3D.controlPts(:,1)-p17(1))<eps),...
                find(abs(skinMesh3D.controlPts(:,2)-p17(2))<eps));


bb1 = intersect(find(abs(stifMesh3D1.controlPts(:,1)-p3(1))<eps),...
                find(abs(stifMesh3D1.controlPts(:,2)-p3(2))<eps));
bb2 = intersect(find(abs(stifMesh3D1.controlPts(:,1)-p5(1))<eps),...
                find(abs(stifMesh3D1.controlPts(:,2)-p5(2))<eps));
bb3 = intersect(find(abs(stifMesh3D1.controlPts(:,1)-p7(1))<eps),...
                find(abs(stifMesh3D1.controlPts(:,2)-p7(2))<eps));            
bb4 = intersect(find(abs(stifMesh3D1.controlPts(:,1)-p9(1))<eps),...
                find(abs(stifMesh3D1.controlPts(:,2)-p9(2))<eps));
            
bb5 = intersect(find(abs(stifMesh3D2.controlPts(:,1)-p11(1))<eps),...
                find(abs(stifMesh3D2.controlPts(:,2)-p11(2))<eps));
bb6 = intersect(find(abs(stifMesh3D2.controlPts(:,1)-p13(1))<eps),...
                find(abs(stifMesh3D2.controlPts(:,2)-p13(2))<eps));
bb7 = intersect(find(abs(stifMesh3D2.controlPts(:,1)-p15(1))<eps),...
                find(abs(stifMesh3D2.controlPts(:,2)-p15(2))<eps));
bb8 = intersect(find(abs(stifMesh3D2.controlPts(:,1)-p17(1))<eps),...
                find(abs(stifMesh3D2.controlPts(:,2)-p17(2))<eps));
            
aa  = [];
bb  = [];

for i=1:length(aa1)
    s = aa1(i);
    e = aa2(i);
    aa = [aa; s:e];
    s = aa3(i);
    e = aa4(i);
    aa = [aa; s:e];
%     s = aa5(i);
%     e = aa6(i);
%     aa = [aa; s:e];
%     s = aa7(i);
%     e = aa8(i);
%     aa = [aa; s:e];
    
    s = bb1(i) + size(skinMesh3D.controlPts,1);
    e = bb2(i) + size(skinMesh3D.controlPts,1);
    bb = [bb; s:e];
    s = bb3(i) + size(skinMesh3D.controlPts,1);
    e = bb4(i) + size(skinMesh3D.controlPts,1);
    bb = [bb; s:e];
%     s = bb5(i) + size(skinMesh3D.controlPts,1) + size(stifMesh3D1.controlPts,1);
%     e = bb6(i) + size(skinMesh3D.controlPts,1) + size(stifMesh3D1.controlPts,1);
%     bb = [bb; s:e];
%     s = bb7(i) + size(skinMesh3D.controlPts,1) + size(stifMesh3D1.controlPts,1);
%     e = bb8(i) + size(skinMesh3D.controlPts,1) + size(stifMesh3D1.controlPts,1);
%     bb = [bb; s:e];
end

aa=aa(:);
bb=bb(:);

bb = bb + size(skinMesh3D.controlPts,1);

for e=1:size(stifMesh3D1.globElems,1)
    sctr  = stifMesh3D1.locElems(e,:);
    sctrg = stifMesh3D1.globElems(e,:);
    [C,IA,IB] = intersect(sctr,bb);
    remainIds = setdiff(1:length(sctr),IA);
    sctr(IA)  = aa(IB);
    sctr(remainIds)  = sctrg(remainIds);
    stifMesh3D1.globElems(e,:) = sctr;
end

% store the meshes in data structure 'mesh'
mesh{1} = skinMesh3D;
mesh{2} = stifMesh3D1;
%mesh{3} = stifMesh3D2;

vmesh{1} = vMeshSkin3D;
vmesh{2} = vMeshStif3D1;
%vmesh{3} = vMeshStif3D2;

%% datastructure, 'data', contains the whole mesh

data.mesh     = mesh;
data.vmesh    = vmesh;
data.pntCount = 0;

for ip=1:length(data.mesh)
  data.pntCount = data.pntCount + size(mesh{ip}.controlPts,1);
end


%% Materials

e1 = 115e3;
e2 = 8.5e3;
e3 = 8.5e3;
nu12 = 0.29;
nu23 = 0.29;
nu31 = 0.29;
g12  = 4.5e3;
g23  = e2/(2+2*nu23);
g31  = 4.5e3;

dirs = [0,90,0,90];

for i =1:length(dirs)
    material = createOrthotropicMaterial (dirs(i),3,...
        e1,e2,e3,nu12,nu23,nu31,g12,g23,g31);
    
    materials{i} = material;
end

%% element sets for each material

for ip=1:length(data.mesh)        % loop over patch
    mesh = data.mesh{ip};         % mesh of the current patch
    noElemsU = mesh.noElemsU;
    noElemsV = mesh.noElemsV;
    noElemsW = mesh.noElemsW;

    xx0 = [1:noElemsU              noElemsU*7+1:noElemsU*noElemsV];
    xx1 = [noElemsU+1  :2*noElemsU noElemsU*6+1:noElemsU*7];
    xx2 = [noElemsU*2+1:noElemsU*3 noElemsU*5+1:noElemsU*6];
    xx3 = [noElemsU*3+1:noElemsU*4 noElemsU*4+1:noElemsU*5];

    cc = xx0;
    for iz=1:noElemsW-1
        cc = [cc xx0 + noElemsU*noElemsV*iz];
    end
    elem{1} = cc;

    cc = xx1;
    for iz=1:noElemsW-1
        cc = [cc xx1 + noElemsU*noElemsV*iz];
    end
    elem{2} = cc;

    cc = xx2;
    for iz=1:noElemsW-1
        cc = [cc xx2 + noElemsU*noElemsV*iz];
    end
    elem{3} = cc;

    cc = xx3;
    for iz=1:noElemsW-1
        cc = [cc xx3 + noElemsU*noElemsV*iz];
    end
    elem{4} = cc;

    elemSet{ip} = elem;
end

tic;

%% find boundary nodes for boundary conditions



%% boundary nodes for traction

noPtsX     = data.mesh{1}.noPtsX;
noPtsY     = data.mesh{1}.noPtsY;
noPtsZ     = data.mesh{1}.noPtsZ;

bndPts1 = zeros(noPtsY,noPtsZ);
bndPts2 = zeros(noPtsY,noPtsZ);

for iy=1:noPtsY
    start1 = noPtsX*iy;
    start2 = noPtsX*(iy-1)+1;
    bndPts1(iy,1) =start1;
    bndPts2(iy,1) =start2;
    for iz=1:noPtsZ-1
        bndPts1(iy,iz+1) =start1 + noPtsX*noPtsY*iz;
        bndPts2(iy,iz+1) =start2 + noPtsX*noPtsY*iz;
    end
end

[bndElems1,index1] = surfaceMesh (data.mesh{1}.wKnot,data.mesh{1}.vKnot,...
    bndPts1,data.mesh{1}.r,data.mesh{1}.q,...
    data.mesh{1}.noPtsZ,data.mesh{1}.noPtsY,...
    data.mesh{1}.rangeW,data.mesh{1}.rangeV,data.mesh{1}.elConnW,data.mesh{1}.elConnV);

[bndElems2,index2] = surfaceMesh (data.mesh{1}.wKnot,data.mesh{1}.vKnot,...
    bndPts2,data.mesh{1}.r,data.mesh{1}.q,...
    data.mesh{1}.noPtsZ,data.mesh{1}.noPtsY,...
    data.mesh{1}.rangeW,data.mesh{1}.rangeV,data.mesh{1}.elConnW,data.mesh{1}.elConnV);

xnodes = unique(bndPts2);
ynodes = xnodes;
znodes = xnodes;

fnodes = unique(bndPts1);

data.xnodes   = xnodes;
data.ynodes   = ynodes;
data.znodes   = znodes;

dirichlet.xnodes = xnodes;
dirichlet.ynodes = xnodes;
dirichlet.znodes = xnodes;

newmann.xnodes   = fnodes;

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



% E0           = 3e6;  % Young modulus
% nu0          = 0.3;  % Poisson’s ratio
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % COMPUTE COMPLIANCE MATRIX
% D=zeros(6,6);
% D(1:3,1:3)=E0/(1+nu0)/(1-2*nu0)*[ 1-nu0 nu0 nu0;
%                                   nu0 1-nu0 nu0;
%                                   nu0 nu0 1-nu0];
% D(4:6,4:6)=E0/2/(1+nu0)*eye(3);
% 
% materials{1}.stiffMat=D;
