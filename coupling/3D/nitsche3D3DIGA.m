% This file implements the Nitsche method to join two 3D meshes.
%
% Discretisation: NURBS finite elements.
%
% Problem: Timoshenko beam in bending.
%
% Vinh Phu Nguyen
% Cardiff University, UK
% 5 July 2013

addpath ../../fem_util/
addpath ../../gmshFiles/
addpath ../../post-processing/
addpath ../../fem-functions/
addpath ../../analytical-solutions/
addpath ../../meshing/
addpath ../../nurbs-util/
addpath ../../nurbs-geopdes/inst/
addpath ../../delamination/
addpath ../../C_files/

clear all
colordef white
state = 0;
tic;

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,'splinecurve.eps',opts)


% MATERIAL PROPERTIES
E0  = 1000;  % Young modulus
nu0 = 0.3;   % Poisson ratio

% BEAM PROPERTIES
L     = 10;     % length of the beam
width = 1;
t     = 1;     % thicknes

% TIP LOAD
P = -1;
q  = 0;    % uniformly distributed loads
ubar=-1;
% Constitutive matrices:
C=zeros(6,6);
C(1:3,1:3)=E0/(1+nu0)/(1-2*nu0)*[ 1-nu0 nu0 nu0;
    nu0 1-nu0 nu0;
    nu0 nu0 1-nu0];
C(4:6,4:6)=E0/2/(1+nu0)*eye(3);

% penalty parameter

alpha = 1e9;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GENERATE FINITE ELEMENT MESH
%
plotMesh  = 1;
disp([num2str(toc),'   GENERATING MESH'])

% MESH PROPERTIES

elemType  = 'B8'; % the element type for solid1;
elemTypeB = 'Q4'; % the element type for surface;

%% domain 1 ----------------------------------------------------------------
L1    = 5;

controlPts          = zeros(4,2,2,2);

controlPts(1:3,1,1,1) = [0;0;0];
controlPts(1:3,2,1,1) = [L1;0;0];
controlPts(1:3,1,2,1) = [0;width;0];
controlPts(1:3,2,2,1) = [L1;width;0];

controlPts(1:3,1,1,2) = [0;0;t];
controlPts(1:3,2,1,2) = [L1;0;t];
controlPts(1:3,1,2,2) = [0;width;t];
controlPts(1:3,2,2,2) = [L1;width;t];

controlPts(4,:,:)   = 1;

uKnot = [0 0 1 1];
vKnot = [0 0 1 1];
wKnot = [0 0 1 1];

solid1 = nrbmak(controlPts,{uKnot vKnot wKnot});

% evaluate order

solid1 = nrbdegelev(solid1,[2 2 2]);


uKnot     = cell2mat(solid1.knots(1));
vKnot     = cell2mat(solid1.knots(2));
wKnot     = cell2mat(solid1.knots(3));

refineCountX = 4;
for i=1:refineCountX
    uKnotVectorU = unique(uKnot);
    uKnotVectorV = unique(vKnot);
    uKnotVectorW = unique(wKnot);
    newKnotsY = [];
    newKnotsX = uKnotVectorU(1:end-1) + 0.5*diff(uKnotVectorU);
    newKnotsY = uKnotVectorV(1:end-1) + 0.5*diff(uKnotVectorV);
    newKnots  = {newKnotsX [] []};
    solid1     = nrbkntins(solid1,newKnots);
    uKnot     = cell2mat(solid1.knots(1));
    vKnot     = cell2mat(solid1.knots(2));
    wKnot     = cell2mat(solid1.knots(3));
end

refineCountX = 2;
for i=1:refineCountX
    uKnotVectorU = unique(uKnot);
    uKnotVectorV = unique(vKnot);
    uKnotVectorW = unique(wKnot);
    
    newKnotsZ = uKnotVectorW(1:end-1) + 0.5*diff(uKnotVectorW);
    newKnotsY = uKnotVectorV(1:end-1) + 0.5*diff(uKnotVectorV);
    newKnots  = {[] newKnotsY newKnotsZ};
    solid1     = nrbkntins(solid1,newKnots);
    uKnot     = cell2mat(solid1.knots(1));
    vKnot     = cell2mat(solid1.knots(2));
    wKnot     = cell2mat(solid1.knots(3));
end

mesh1      = buildIGA3DMesh (solid1);

node1     = mesh1.controlPts;
element1  = mesh1.locElems;

numx1 = mesh1.noElemsU;
numy1 = mesh1.noElemsV;
numz1 = mesh1.noElemsW;

% domain 2 ----------------------------------------------------------------

controlPts          = zeros(4,2,2,2);

controlPts(1:3,1,1,1) = [L1;0;0];
controlPts(1:3,2,1,1) = [L;0;0];
controlPts(1:3,1,2,1) = [L1;width;0];
controlPts(1:3,2,2,1) = [L;width;0];

controlPts(1:3,1,1,2) = [L1;0;t];
controlPts(1:3,2,1,2) = [L;0;t];
controlPts(1:3,1,2,2) = [L1;width;t];
controlPts(1:3,2,2,2) = [L;width;t];

controlPts(4,:,:)   = 1;

uKnot = [0 0 1 1];
vKnot = [0 0 1 1];
wKnot = [0 0 1 1];

solid2 = nrbmak(controlPts,{uKnot vKnot wKnot});

% evaluate order

solid2 = nrbdegelev(solid2,[2 2 2]);


uKnot     = cell2mat(solid2.knots(1));
vKnot     = cell2mat(solid2.knots(2));

refineCountX = 4;
for i=1:refineCountX
    uKnotVectorU = unique(uKnot);
    uKnotVectorV = unique(vKnot);
    newKnotsY = [];
    newKnotsX = uKnotVectorU(1:end-1) + 0.5*diff(uKnotVectorU);
    newKnotsY = uKnotVectorV(1:end-1) + 0.5*diff(uKnotVectorV);
    newKnots  = {newKnotsX [] []};
    solid2    = nrbkntins(solid2,newKnots);
    uKnot     = cell2mat(solid2.knots(1));
    vKnot     = cell2mat(solid2.knots(2));
end

newKnots  = {[] [] [0.5]};
solid2    = nrbkntins(solid2,newKnots);

mesh2     = buildIGA3DMesh (solid2);

node2     = mesh2.controlPts;
element2  = mesh2.locElems;
mesh2.globElems=mesh2.globElems+size(node1,1);

numx2 = mesh2.noElemsU;
numy2 = mesh2.noElemsV;
numz2 = mesh2.noElemsW;

%% ------------------------------------------------------------------------
% boundary mesh where coupling terms are evaluated

bndNodes  = find(abs(node1(:,1)-L1)<1e-14);

[bndMesh1,botIndex]     = surfaceMesh (mesh1.vKnot,mesh1.wKnot,bndNodes,mesh1.q,mesh1.r,...
    mesh1.noPtsY,mesh1.noPtsZ,mesh1.rangeV,mesh1.rangeW,mesh1.elConnV,mesh1.elConnW);

%map1 = [numx1 numx1*2 numx1*3 numx1*4];

map1 = [];

for iz=1:numz1
    for iy=1:numy1
        map1 = [map1;numx1*iy + numx1*numy1*(iz-1)];
    end
end
% 
map2 = [1 1 1 1 1 1 1 1 ...
    numx2+1 numx2+1 numx2+1 numx2+1 numx2+1 numx2+1 numx2+1 numx2+1];

%map2=1;

ngp = (mesh1.q+1) * (mesh1.r+1);

[W1,Q1]=quadrature( mesh1.q+1, 'GAUSS', 2 ); % two point quadrature
%W1=1;Q1=[0 0];

GP1 = [];
GP2 = [];

for e=1:size(bndMesh1,1)
    be1      = map1(e);
    be2      = map2(e);
    sctrS    = bndMesh1(e,:);    % connectivity of surface
    sctrV1   = element1(be1,:);  % connectivity of solid1
    sctrV2   = element2(be2,:);  % connectivity of solid2
    pts1     = node1(sctrV1,:);  % node coords of element1
    pts2     = node2(sctrV2,:);  % node coords of element2
    ptsS     = node1(sctrS,:);   % node coords of coupling surface
    
    idv    = botIndex(e,1);
    idw    = botIndex(e,2);
    
    etaE   = mesh1.rangeV(idv,:); % [eta_j,eta_j+1]
    zetaE  = mesh1.rangeW(idw,:); % [zeta_k,zeta_k+1]
    
    xiE1     = mesh1.rangeU(end,:);                  % [xi_i,xi_i+1]
    etE1     = mesh1.rangeV(mesh1.index(be1,2),:);   % [et_j,eta_j+1]
    zetE1    = mesh1.rangeW(mesh1.index(be1,3),:);   % [zet_j,zeta_j+1]
    
    xiE2     = mesh2.rangeU(1,:);                    % [xi_i,xi_i+1]
    etE2     = mesh2.rangeV(mesh2.index(be2,2),:);   % [et_j,eta_j+1]
    zetE2    = mesh2.rangeW(mesh2.index(be2,3),:);   % [zet_j,zeta_j+1]
    
    for q=1:size(W1)
        pt=Q1(q,:);
        wt=W1(q);
        
        Xi      = parent2ParametricSpace(etaE, pt(1));
        Eta     = parent2ParametricSpace(zetaE,pt(2));
        J2      = jacobianPaPaMapping(etaE,zetaE);
        
        % compute derivatives of shape functions
        [R dRdxi dRdeta] = NURBS2DBasisDers([Xi; Eta],mesh1.q,mesh1.r,mesh1.vKnot,mesh1.wKnot,mesh1.weights');
        J0       = [dRdxi; dRdeta]*ptsS;
        a1       = J0(1,:);
        a2       = J0(2,:);
        a3       = cross(a1,a2);
        norma    = norm(a3); a3 = a3/norma;
        
        x   = R*ptsS;
        X1  = global2LocalMapNURBS3D(x,pts1,xiE1,etE1,zetE1,mesh1.p,mesh1.q,mesh1.r,...
            mesh1.uKnot,mesh1.vKnot,mesh1.wKnot,mesh1.weights);
        X2  = global2LocalMapNURBS3D(x,pts2,xiE2,etE2,zetE2,mesh2.p,mesh2.q,mesh2.r,...
            mesh2.uKnot,mesh2.vKnot,mesh2.wKnot,mesh2.weights);
        GP1 = [GP1;X1 wt*norma*J2 a3];
        GP2 = [GP2;X2];
    end
end

%GP1(:,1) = 1;

%% ------------------------------------------------------------------------
% boundary nodes

fixedNode = find(node1(:,1)==0);
rightNode = find(node2(:,1)==L);
rightNode = rightNode + size(node1,1);

uNode = [fixedNode; rightNode];

uFixed    = zeros(1,length(fixedNode));  % a vector of u_x for the nodes
vFixed    = zeros(1,length(fixedNode));
wFixed    = [zeros(1,length(fixedNode)) ubar*ones(1,length(rightNode))];

udofs     = 3*fixedNode-2;
vdofs     = 3*fixedNode-1;
wdofs     = 3*uNode-0;

numnodes = size(node1,1) + size(node2,1);

%% find active dofs used in eigenvalue problem
activeNodes=setdiff(1:numnodes,uNode);

activeDofs=[];
for i=1:length(activeNodes)
    in = activeNodes(i);
    activeDofs=[activeDofs; 3*in-2;3*in-1;3*in];
end

% PLOT MESH
if ( plotMesh )
    figure
    hold on
    nrbkntplot(solid1);
    nrbkntplot(solid2);
    plot3(mesh1.controlPts(:,1),mesh1.controlPts(:,2),mesh1.controlPts(:,3),'o','MarkerEdgeColor','k',...
        'MarkerFaceColor','r','MarkerSize',4.5);
    plot3(mesh2.controlPts(:,1),mesh2.controlPts(:,2),mesh2.controlPts(:,3),'s','MarkerEdgeColor','k',...
        'MarkerFaceColor','r','MarkerSize',4.5);
    axis off
    %plot3(node1(fixedNode,1),node1(fixedNode,2),node1(fixedNode,3),'r*');
    %    plot3(node2(rightNode,1),node2(rightNode,2),node2(rightNode,3),'rs');
    axis off
    %axis([0 L -c c])
end

%% Combine meshes into one data structure

matMap1 = ones(mesh1.noElems,1);
matMap2 = ones(mesh2.noElems,1);

vMesh1=buildVisualizationMesh3D(solid1);
vMesh2=buildVisualizationMesh3D(solid2);

meshes{1}  = mesh1;
meshes{2}  = mesh2;
vmeshes{1} = vMesh1;
vmeshes{2} = vMesh2;


meshData.mesh    = meshes;
meshData.vmesh   = vmeshes;
% meshData.noElems = meshes{1}.noElems;
% meshData.noPts   = meshes{1}.noPts;
meshData.matMap{1}=matMap1;
meshData.matMap{2}=matMap2;
material.stiffMat=C;
materials{1}      = material;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE SYSTEM DATA STRUCTURES

f=zeros(3*numnodes,1);          % external load vector
K=zeros(3*numnodes,3*numnodes);  % stiffness matrix

% ******************************************************************************
% ***                          P R O C E S S I N G                           ***
% ******************************************************************************
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%% COMPUTE STIFFNESS MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp([num2str(toc),'   COMPUTING STIFFNESS MATRIX'])

[W,Q]=quadrature( mesh1.p+1, 'GAUSS', 3 ); % 2x2x2 Gaussian quadrature

%% for domain 1 (continuum)------------------------------------------------

for ip=1:length(meshes)                               % start of patch loop
    mesh = meshes{ip};
    elementg = mesh.globElems;
    elementl = mesh.locElems;
    node     = mesh.controlPts;
    for e=1:size(elementg,1)                          % start of element loop
        sctrl           = elementl(e,:);              % element scatter vector
        sctrg           = elementg(e,:);              % element scatter vector
        nn              = length(sctrg);
        sctrB(1:3:3*nn) = 3*sctrg-2;
        sctrB(2:3:3*nn) = 3*sctrg-1;
        sctrB(3:3:3*nn) = 3*sctrg-0;
        pts             = node(sctrl,:);
        
        idu    = mesh.index(e,1);
        idv    = mesh.index(e,2);
        idw    = mesh.index(e,3);
        
        xiE    = mesh.rangeU(idu,:); % [xi_i,xi_i+1]
        etaE   = mesh.rangeV(idv,:); % [eta_j,eta_j+1]
        zetaE  = mesh.rangeW(idw,:); % [zeta_k,zeta_k+1]
        
        for q=1:size(W,1)
            pt=Q(q,:);
            wt=W(q);
            
            Xi      = parent2ParametricSpace(xiE,  pt(1));
            Eta     = parent2ParametricSpace(etaE, pt(2));
            Zeta    = parent2ParametricSpace(zetaE,pt(3));
            J2      = jacobianPaPaMapping3d (xiE,etaE,zetaE);
            % compute derivative of basis functions w.r.t parameter coord
            [N dRdxi dRdeta dRdzeta] = NURBS3DBasisDers([Xi;Eta;Zeta],...
                mesh.p,mesh.q,mesh.r,mesh.uKnot,mesh.vKnot,mesh.wKnot,mesh.weights');
            jacob      = pts'*[dRdxi' dRdeta' dRdzeta'];
            J1         = det(jacob);
            invJacob   = inv(jacob);
            dRdx       = [dRdxi' dRdeta' dRdzeta'] * invJacob;
            B          = getBmatrix3D(nn,dRdx);
            K(sctrB,sctrB) = K(sctrB,sctrB) + B' * C * B * J1 * J2 * wt;
        end  % of quadrature loop
    end
end

%% determination of alpha--------------------------------------------------

Ktilde = K;
H=zeros(size(K,1),size(K,2)); % stiffness matrix

for i=1:size(bndMesh1,1)                     % start of element loop
    e1     = map1(i);
    e2     = map2(i);
    sctr1  = element1(e1,:);
    sctr2  = element2(e2,:);
    sctr2n = sctr2 + size(node1,1);
    
    nn1    = length(sctr1);
    nn2    = length(sctr2);
    
    sctrB1(1:3:3*nn1) = 3*sctr1-2;
    sctrB1(2:3:3*nn1) = 3*sctr1-1;
    sctrB1(3:3:3*nn1) = 3*sctr1-0;
    
    sctrB2(1:3:3*nn2) = 3*sctr2n-2;
    sctrB2(2:3:3*nn2) = 3*sctr2n-1;
    sctrB2(3:3:3*nn2) = 3*sctr2n-0;
    
    pts1 = node1(sctr1,:);
    pts2 = node2(sctr2,:);
    
    Kp11 = zeros(nn1*3,nn1*3);
    Kp12 = zeros(nn1*3,nn2*3);
    Kp22 = zeros(nn2*3,nn2*3);
    
    Kd11 = zeros(nn1*3,nn1*3);
    Kd12 = zeros(nn1*3,nn2*3);
    Kd21 = zeros(nn2*3,nn1*3);
    Kd22 = zeros(nn2*3,nn2*3);
    
    for q=ngp*(i-1)+1:ngp*(i-1)+ngp
        pt1    = GP1(q,1:3);
        pt2    = GP2(q,1:3);
        wt1    = GP1(q,4);
        
        normal = GP1(q,5:end);
        n = [normal(1) 0 0 normal(2) 0 normal(3);...
            0 normal(2) 0 normal(1) normal(3) 0;...
            0 0 normal(3) 0 normal(2) normal(1)];
        %n=-n;
        
        [N1 dN1dxi dN1deta dN1dzeta] = NURBS3DBasisDers(pt1,...
            mesh1.p,mesh1.q,mesh1.r,mesh1.uKnot,mesh1.vKnot,mesh1.wKnot,mesh1.weights');
        [N2 dN2dxi dN2deta dN2dzeta] = NURBS3DBasisDers(pt2,...
            mesh2.p,mesh2.q,mesh2.r,mesh2.uKnot,mesh2.vKnot,mesh2.wKnot,mesh2.weights');
        
        J1      = pts1'*[dN1dxi' dN1deta' dN1dzeta'];
        J2      = pts2'*[dN2dxi' dN2deta' dN2dzeta'];
        
        dN1dx   = [dN1dxi' dN1deta' dN1dzeta']*inv(J1);
        dN2dx   = [dN2dxi' dN2deta' dN2dzeta']*inv(J2);
        
        B1      = getBmatrix3D(nn1,dN1dx);
        B2      = getBmatrix3D(nn2,dN2dx);
                                  
        Kd11 = Kd11 + B1'*C'*n'*n*C*B1 *wt1;
        Kd12 = Kd12 + B1'*C'*n'*n*C*B2 *wt1;
        Kd21 = Kd21 + B2'*C'*n'*n*C*B1 *wt1;
        Kd22 = Kd22 + B2'*C'*n'*n*C*B2 *wt1;   
        
    end  % of quadrature loop
    
    H(sctrB1,sctrB1)  = H(sctrB1,sctrB1)  + Kd11;
    H(sctrB1,sctrB2)  = H(sctrB1,sctrB2)  + Kd12;
    H(sctrB2,sctrB1)  = H(sctrB2,sctrB1)  + Kd21;
    H(sctrB2,sctrB2)  = H(sctrB2,sctrB2)  + Kd22;
end

lambda = eigs(inv(Ktilde(activeDofs,activeDofs))*H(activeDofs,activeDofs),5,'lm');

alpha=max(lambda)/2

%% interface integrals-----------------------------------------------------

for i=1:size(bndMesh1,1)                     % start of element loop
    e1     = map1(i);
    e2     = map2(i);
    sctr1  = element1(e1,:);
    sctr2  = element2(e2,:);
    sctr2n = sctr2 + size(node1,1);
    
    nn1    = length(sctr1);
    nn2    = length(sctr2);
    
    sctrB1(1:3:3*nn1) = 3*sctr1-2;
    sctrB1(2:3:3*nn1) = 3*sctr1-1;
    sctrB1(3:3:3*nn1) = 3*sctr1-0;
    
    sctrB2(1:3:3*nn2) = 3*sctr2n-2;
    sctrB2(2:3:3*nn2) = 3*sctr2n-1;
    sctrB2(3:3:3*nn2) = 3*sctr2n-0;
    
    pts1 = node1(sctr1,:);
    pts2 = node2(sctr2,:);
    
    Kp11 = zeros(nn1*3,nn1*3);
    Kp12 = zeros(nn1*3,nn2*3);
    Kp22 = zeros(nn2*3,nn2*3);
    
    Kd11 = zeros(nn1*3,nn1*3);
    Kd12 = zeros(nn1*3,nn2*3);
    Kd21 = zeros(nn2*3,nn1*3);
    Kd22 = zeros(nn2*3,nn2*3);
    
    for q=ngp*(i-1)+1:ngp*(i-1)+ngp
        pt1    = GP1(q,1:3);
        pt2    = GP2(q,1:3);
        wt1    = GP1(q,4);
        
        normal = GP1(q,5:end);
        n = [normal(1) 0 0 normal(2) 0 normal(3);...
            0 normal(2) 0 normal(1) normal(3) 0;...
            0 0 normal(3) 0 normal(2) normal(1)];
        %n=-n;
        
        [N1 dN1dxi dN1deta dN1dzeta] = NURBS3DBasisDers(pt1,...
            mesh1.p,mesh1.q,mesh1.r,mesh1.uKnot,mesh1.vKnot,mesh1.wKnot,mesh1.weights');
        [N2 dN2dxi dN2deta dN2dzeta] = NURBS3DBasisDers(pt2,...
            mesh2.p,mesh2.q,mesh2.r,mesh2.uKnot,mesh2.vKnot,mesh2.wKnot,mesh2.weights');
        
        J1      = pts1'*[dN1dxi' dN1deta' dN1dzeta'];
        J2      = pts2'*[dN2dxi' dN2deta' dN2dzeta'];
        
        dN1dx   = [dN1dxi' dN1deta' dN1dzeta']*inv(J1);
        dN2dx   = [dN2dxi' dN2deta' dN2dzeta']*inv(J2);
        
        B1      = getBmatrix3D(nn1,dN1dx);
        B2      = getBmatrix3D(nn2,dN2dx);
        
        Nm1(1,1:3:3*nn1)  = N1';
        Nm1(2,2:3:3*nn1)  = N1';
        Nm1(3,3:3:3*nn1)  = N1';
        
        Nm2(1,1:3:3*nn2)  = N2';
        Nm2(2,2:3:3*nn2)  = N2';
        Nm2(3,3:3:3*nn2)  = N2';
        
        % COMPUTE ELEMENT STIFFNESS AT QUADRATURE POINT
        
        Kp11 = Kp11 + alpha*(Nm1'*Nm1)*wt1;
        Kp12 = Kp12 + alpha*(Nm1'*Nm2)*wt1;
        Kp22 = Kp22 + alpha*(Nm2'*Nm2)*wt1;
        
        Kd11 = Kd11 + 0.5 * Nm1'* n * C * B1 *wt1;
        Kd12 = Kd12 + 0.5 * Nm1'* n * C * B2 *wt1;
        Kd21 = Kd21 + 0.5 * Nm2'* n * C * B1 *wt1;
        Kd22 = Kd22 + 0.5 * Nm2'* n * C * B2 *wt1;
        
    end  % of quadrature loop
    
    K(sctrB1,sctrB1)  = K(sctrB1,sctrB1)  - Kd11 - Kd11' + Kp11;
    K(sctrB1,sctrB2)  = K(sctrB1,sctrB2)  - Kd12 + Kd21' - Kp12;
    K(sctrB2,sctrB1)  = K(sctrB2,sctrB1)  + Kd21 - Kd12' - Kp12';
    K(sctrB2,sctrB2)  = K(sctrB2,sctrB2)  + Kd22 + Kd22' + Kp22;
end


%f(3*rightNode) = P;

%%%%%%%%%%%%%%%%%%% END OF STIFFNESS MATRIX COMPUTATION %%%%%%%%%%%%%%%%%%%


% APPLY ESSENTIAL BOUNDARY CONDITIONS
disp([num2str(toc),'   APPLYING BOUNDARY CONDITIONS'])

[K,f]=applyDirichletBCs3D(K,f,udofs,vdofs,wdofs,uFixed,vFixed,wFixed);

% SOLVE SYSTEM
disp([num2str(toc),'   SOLVING SYSTEM'])
U=K\f;

%******************************************************************************
%*** POST - PROCESSING ***
%***************************************************

numnode1 = size(node1,1);
numnode2 = size(node2,1);

index1   = 1:numnode1;
index2   = numnode1+1:numnodes;

Ux1      = U(3*index1-2);
Uy1      = U(3*index1-1);
Uz1      = U(3*index1-0);

Ux2      = U(3*index2-2);
Uy2      = U(3*index2-1);
Uz2      = U(3*index2-0);

% Here we plot the stresses and displacements of the solution.
disp([num2str(toc),'   POST-PROCESSING'])

dispNorm=L/max(sqrt(Ux1.^2+Uy1.^2));
scaleFact=0.05*dispNorm;

element1(:,[3 4 7 8]) = element1(:,[4 3 8 7]);
element2(:,[3 4 7 8]) = element2(:,[4 3 8 7]);
colordef white
figure
clf
hold on
%plot_field(node1+scaleFact*[Ux1 Uy1 Uz1],element1,elemType,Uy1);
%plot_field(node2+scaleFact*[Ux2 Uy2 Uz2],element2,elemType,Uy2);
plot_mesh(node1+scaleFact*[Ux1 Uy1 Uz1],element1,elemType,'r-',1);
plot_mesh(node2+scaleFact*[Ux2 Uy2 Uz2],element2,elemType,'b-',1);
%colorbar
axis off
%title('DEFORMED DISPLACEMENT IN Y-DIRECTION')

%%

displacement = [U(1:3:end) U(2:3:end) U(3:3:end)];
damage       = zeros(length(U),1);


vtuFileName = 'nitsche-3D3D-iga';
for ip=1:2
    vtuFile = strcat(vtuFileName,'-mesh',num2str(ip));
    %figure; hold on;
    ok      = plotStress3DForPatch(meshData,ip,vtuFile,displacement,damage,materials);
end

pvdFile = fopen(strcat(vtuFileName,'.pvd'), 'wt');

fprintf(pvdFile,'<VTKFile byte_order="LittleEndian" type="Collection" version="0.1">\n');
fprintf(pvdFile,'<Collection>\n');


for im=1:length(meshes)
    vtuFile = sprintf('%s%s%d%s',vtuFileName,'-mesh',im,'.vtu');
    fprintf(pvdFile,'<DataSet file=''%s'' groups='''' part=''%d'' timestep=''%d''/>\n',...
        vtuFile,im,1);
end


fprintf(pvdFile,'</Collection>\n');
fprintf(pvdFile,'</VTKFile>\n');

fclose(pvdFile);


