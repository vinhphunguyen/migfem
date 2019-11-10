% This file implements the Nitsche method to join two mechanical models:
% a 3D continuum model and a Kirchhoff plate model.
%
% Discretisation: NURBS elements.
%
% Non-conforming coupling.
% A square Kirchhoff plate is overlapped in the middle by a 3D solid.
% This implementation illustrates the case of multiple (4) coupling
% surfaces. 
%
% Vinh Phu Nguyen
% Cardiff University, UK
% 20 July 2013

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
addpath ../

clear all
colordef white
state = 0;
tic;

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,'splinecurve.eps',opts)


% MATERIAL PROPERTIES
E0  = 1000;  % Young modulus
nu0 = 0.3;   % Poisson ratio
k   = 5/6;   % shear correction factor

% BEAM PROPERTIES
L     = 400;     % length of the beam
width = 400;     % width
t     = 20;      % thicknes
L1    = 85;      % size of the solid model

q0    = -10;    % uniformly distributed loads

% Constitutive matrices: continuum and plate
C=zeros(6,6);
C(1:3,1:3)=E0/(1+nu0)/(1-2*nu0)*[ 1-nu0 nu0 nu0;
    nu0 1-nu0 nu0;
    nu0 nu0 1-nu0];
C(4:6,4:6)=E0/2/(1+nu0)*eye(3);

D   = E0*t^3/(12*(1-nu0^2));
Cp  = D*[1 nu0 0;nu0 1 0; 0 0 0.5*(1-nu0)];

% penalty parameter

alpha = 1e6;

EPS = 1e-10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GENERATE FINITE ELEMENT MESH
%
plotMesh  = 1;
disp([num2str(toc),'   GENERATING MESH'])

% MESH PROPERTIES

%% domain 1 ---------------------------------------------------------------

controlPts          = zeros(4,2,2,2);

controlPts(1:3,1,1,1) = [0;0;0];
controlPts(1:3,2,1,1) = [L1;0;0];
controlPts(1:3,1,2,1) = [0;L1;0];
controlPts(1:3,2,2,1) = [L1;L1;0];

controlPts(1:3,1,1,2) = [0;0;t];
controlPts(1:3,2,1,2) = [L1;0;t];
controlPts(1:3,1,2,2) = [0;L1;t];
controlPts(1:3,2,2,2) = [L1;L1;t];

controlPts(4,:,:)   = 1;

uKnot = [0 0 1 1];
vKnot = [0 0 1 1];
wKnot = [0 0 1 1];

solid1 = nrbmak(controlPts,{uKnot vKnot wKnot});

% evaluate order

solid1 = nrbdegelev(solid1,[2 2 2]);

uKnot     = cell2mat(solid1.knots(1));
vKnot     = cell2mat(solid1.knots(2));

refineCountX = 3;
for i=1:refineCountX
    uKnotVectorU = unique(uKnot);
    uKnotVectorV = unique(vKnot);
    newKnotsY = [];
    newKnotsX = uKnotVectorU(1:end-1) + 0.5*diff(uKnotVectorU);
    newKnotsY = uKnotVectorV(1:end-1) + 0.5*diff(uKnotVectorV);
    newKnots  = {newKnotsX newKnotsY []};
    solid1     = nrbkntins(solid1,newKnots);
    uKnot      = cell2mat(solid1.knots(1));
    vKnot      = cell2mat(solid1.knots(2));
end

solid1     = nrbkntins(solid1,{[] [] [0.5]});

% move the solid to the center of the plate
dx=50;
ctrpts1 = reshape(solid1.coefs,4,[]).';
ctrpts1(:,[1 2 3]) = ctrpts1(:,[1 2 3]) +...
    repmat([0.5*(L-L1)+dx 0.5*(L-L1) -0.5*t],prod(solid1.number),1);
ctrpts1 = reshape(ctrpts1.', 4, solid1.number(1), solid1.number(2),solid1.number(3));
solid1 = nrbmak(ctrpts1,{solid1.knots{1} solid1.knots{2} solid1.knots{3}});

mesh1     = buildIGA3DMesh (solid1);

node1     = mesh1.controlPts;
element1  = mesh1.locElems;

numx1 = mesh1.noElemsU;
numy1 = mesh1.noElemsV;
numz1 = mesh1.noElemsW;

%% domain 2 (plate)--------------------------------------------------------

controlPts          = zeros(4,2,2);

controlPts(1:3,1,1) = [0;0;0];
controlPts(1:3,2,1) = [L;0;0];
controlPts(1:3,1,2) = [0;L;0];
controlPts(1:3,2,2) = [L;L;0];

controlPts(4,:,:)   = 1;

uKnot = [0 0 1 1];
vKnot = [0 0 1 1];

solid2 = nrbmak(controlPts,{uKnot vKnot});

% evaluate order

solid2 = nrbdegelev(solid2,[2 2]);

uKnot     = cell2mat(solid2.knots(1));
vKnot     = cell2mat(solid2.knots(2));

solid2    = nrbkntins(solid2,{[] [0.5] });
refineCountX = 4;
for i=1:refineCountX
    uKnotVectorU = unique(uKnot);
    uKnotVectorV = unique(vKnot);
    newKnotsX = uKnotVectorU(1:end-1) + 0.5*diff(uKnotVectorU);
    newKnotsY = uKnotVectorV(1:end-1) + 0.5*diff(uKnotVectorV);
    newKnots  = {newKnotsX newKnotsY};
    solid2     = nrbkntins(solid2,newKnots);
    uKnot      = cell2mat(solid2.knots(1));
    vKnot      = cell2mat(solid2.knots(2));
end

mesh2     = buildIGA2DMeshForSurface (solid2);

numx2 = mesh2.noElemsU;
numy2 = mesh2.noElemsV;

node2     = mesh2.controlPts;
element2  = mesh2.locElems;
mesh2.globElems=mesh2.globElems+3*size(node1,1);

clamped=1;


%% boundary surfaces-------------------------------------------------------

bndNodesRight  = find(abs(node1(:,1)-0.5*(L+L1)-dx)<EPS);
bndNodesLeft   = find(abs(node1(:,1)-0.5*(L-L1)-dx)<EPS);
bndNodesBot    = find(abs(node1(:,2)-0.5*(L-L1))<EPS);
bndNodesTop    = find(abs(node1(:,2)-0.5*(L+L1))<EPS);

[bndMesh1Right,rightIndex] = surfaceMesh (mesh1.vKnot,mesh1.wKnot,bndNodesRight,mesh1.q,mesh1.r,...
    mesh1.noPtsY,mesh1.noPtsZ,mesh1.rangeV,mesh1.rangeW,mesh1.elConnV,mesh1.elConnW);

[bndMesh1Left,leftIndex] = surfaceMesh (mesh1.vKnot,mesh1.wKnot,bndNodesLeft,mesh1.q,mesh1.r,...
    mesh1.noPtsY,mesh1.noPtsZ,mesh1.rangeV,mesh1.rangeW,mesh1.elConnV,mesh1.elConnW);

[bndMesh1Bot,botIndex] = surfaceMesh (mesh1.uKnot,mesh1.wKnot,bndNodesBot,mesh1.p,mesh1.r,...
    mesh1.noPtsX,mesh1.noPtsZ,mesh1.rangeU,mesh1.rangeW,mesh1.elConnU,mesh1.elConnW);

[bndMesh1Top,topIndex] = surfaceMesh (mesh1.uKnot,mesh1.wKnot,bndNodesTop,mesh1.p,mesh1.r,...
    mesh1.noPtsX,mesh1.noPtsZ,mesh1.rangeU,mesh1.rangeW,mesh1.elConnU,mesh1.elConnW);

if isempty(rightIndex) || isempty(leftIndex) || isempty(botIndex) || isempty(topIndex)
    error('wrong')
end

%% index of surface elements
mapRight = [];
mapLeft  = [];
mapBot   = [];
mapTop   = [];

for iz=1:numz1
    for iy=1:numy1
        mapRight = [mapRight;numx1*iy + numx1*numy1*(iz-1)];
    end
end

for iz=1:numz1
    for iy=1:numy1
        mapLeft = [mapLeft;1+numx1*(iy-1) + numx1*numy1*(iz-1)];
    end
end

for iz=1:numz1
    for iy=1:numx1
        mapBot = [mapBot;iy + numx1*(iz-1)];
    end
end

for iz=1:numz1
    for iy=1:numx1
        mapTop = [mapTop;numx1*(numy1-1)+iy + numx1*numy1*(iz-1)];
    end
end

%% data structure to store the boundary solid elements where they are coupled
%  to plate elements

couplingSurfaces{1}.surfMesh = bndMesh1Right;
couplingSurfaces{1}.elemIDs  = mapRight;
couplingSurfaces{1}.index    = rightIndex;
couplingSurfaces{1}.order    = [mesh1.q mesh1.r];
couplingSurfaces{1}.range{1} = mesh1.rangeV;
couplingSurfaces{1}.range{2} = mesh1.rangeW;
couplingSurfaces{1}.knots{1} = mesh1.vKnot;
couplingSurfaces{1}.knots{2} = mesh1.wKnot;

couplingSurfaces{2}.surfMesh = bndMesh1Left;
couplingSurfaces{2}.elemIDs  = mapLeft;
couplingSurfaces{2}.index    = leftIndex;
couplingSurfaces{2}.order    = [mesh1.q mesh1.r];
couplingSurfaces{2}.range{1} = mesh1.rangeV;
couplingSurfaces{2}.range{2} = mesh1.rangeW;
couplingSurfaces{2}.knots{1} = mesh1.vKnot;
couplingSurfaces{2}.knots{2} = mesh1.wKnot;

couplingSurfaces{3}.surfMesh = bndMesh1Bot;
couplingSurfaces{3}.elemIDs  = mapBot;
couplingSurfaces{3}.index    = botIndex;
couplingSurfaces{3}.order    = [mesh1.p mesh1.r];
couplingSurfaces{3}.range{1} = mesh1.rangeU;
couplingSurfaces{3}.range{2} = mesh1.rangeW;
couplingSurfaces{3}.knots{1} = mesh1.uKnot;
couplingSurfaces{3}.knots{2} = mesh1.wKnot;

couplingSurfaces{4}.surfMesh = bndMesh1Top;
couplingSurfaces{4}.elemIDs  = mapTop;
couplingSurfaces{4}.index    = topIndex;
couplingSurfaces{4}.order    = [mesh1.p mesh1.r];
couplingSurfaces{4}.range{1} = mesh1.rangeU;
couplingSurfaces{4}.range{2} = mesh1.rangeW;
couplingSurfaces{4}.knots{1} = mesh1.uKnot;
couplingSurfaces{4}.knots{2} = mesh1.wKnot;


%% determination of void, cut and standard elements
% using level set representation of  the coupling surface
% which is a square in this example.

% center and dimensions of the intersection surface of solid and plate
x0 = L/2+dx;
y0 = L/2;
a  = L1;
b  = L1;

vMesh=buildVisualizationMesh2D(solid2);

nNode = size(vMesh.node,1);
chi   = zeros(1,nNode);

for i = 1 : nNode
    x  = vMesh.node(i,1);
    y  = vMesh.node(i,2);
    phi1 = -((x-x0)+a/2);
    phi2 =   (x-x0)-a/2;
    phi3 = -((y-y0)+b/2);
    phi4 =   (y-y0)-b/2;
    phix = max(phi1,phi2);
    phiy = max(phi3,phi4);
    chi(1,i)  = max(phix,phiy);
end

elem2State = zeros(size(element2,1),1);

% loop over elements

for iel = 1 : size(element2,1)
    sctr    = element2(iel,:);
    sctrV   = vMesh.element(iel,:);
    
    phi  = chi(1,sctrV);
    
    if     ( max(phi)*min(phi) < 0 )
        elem2State(iel) = 1;             % cut elements
    elseif ( max(phi) < 0 )
        elem2State(iel) = 0;             % void elements
    else
        elem2State(iel) = 2;             % normal elements
    end
end

% find inactive plate dofs

inactiveNodes1 = [];
inactiveNodes2 = [];

for i=1:size(element2,1)
    sctr = element2(i,:);
    if     elem2State(i) == 0
        inactiveNodes1 = [inactiveNodes1;sctr];
    elseif elem2State(i) == 1
        inactiveNodes2 = [inactiveNodes2;sctr];
    end
end

inactiveNodes1 = unique(inactiveNodes1);
inactiveNodes2 = unique(inactiveNodes2);

inactivePlateNodes = setdiff(inactiveNodes1,inactiveNodes2);

voidElems = find(elem2State==0);
cutElems  = find(elem2State==1);

%% loop over coupling surfaces, then solid elements in each surface
% defined GPs and interacting plate elements...

for is=1:length(couplingSurfaces)
    surf    = couplingSurfaces{is};
    sMesh   = surf.surfMesh;
    elemIds = surf.elemIDs;
    index   = surf.index;
    order   = surf.order;
    range1  = surf.range{1};
    range2  = surf.range{2};
    uknot   = surf.knots{1};
    vknot   = surf.knots{2};
       
    [W1,Q1] = quadrature( order(1)+1, 'GAUSS', 2 ); % two point quadrature
    
    for e=1:length(elemIds)
        be       = elemIds (e);
        sctrS    = sMesh(e,:);
        sctrV    = element1(be,:);
        pts1     = node1(sctrV,:);
        ptsS     = node1(sctrS,:);
        
        idv    = index(e,1);
        idw    = index(e,2);
        % parametric range of surface
        etaE   = range1(idv,:); % [eta_j,eta_j+1]
        zetaE  = range2(idw,:); % [zeta_k,zeta_k+1]
        % parametric range of volume
        xiE1   = mesh1.rangeU(mesh1.index(be,1),:);   % [xi_i,xi_i+1]
        etE1   = mesh1.rangeV(mesh1.index(be,2),:);   % [et_j,eta_j+1]
        zetE1  = mesh1.rangeW(mesh1.index(be,3),:);   % [zet_j,zeta_j+1]
        
        GP1 = [];
        GP2 = [];
    
        for q=1:size(W1)
            pt=Q1(q,:);
            wt=W1(q);            
            Xi      = parent2ParametricSpace(etaE, pt(1));
            Eta     = parent2ParametricSpace(zetaE,pt(2));
            J2      = jacobianPaPaMapping(etaE,zetaE);            
            % compute derivatives of shape functions
            [N2, dN2dxi, dN2deta] = NURBS2DBasisDers([Xi;Eta],order(1),order(2),uknot,vknot,mesh1.weights');
            J0       = [dN2dxi; dN2deta]*ptsS;
            a1       = J0(1,:);
            a2       = J0(2,:);
            a3       = cross(a1,a2);  norma    = norm(a3); a3 = a3/norma;
            
            x   = N2*ptsS;
            xP  = x([1 2]);  % project of x on the plate
            X1  = global2LocalMapNURBS3D(x,pts1,xiE1,etE1,zetE1,mesh1.p,mesh1.q,mesh1.r,...
                mesh1.uKnot,mesh1.vKnot,mesh1.wKnot,mesh1.weights);
            
            GP1 = [GP1;X1 wt*norma*J2 a3];
            for ce=1:length(cutElems)
                id   = cutElems(ce);
                sctrV = vMesh.element(id,:);
                sctr = mesh2.locElems(id,:);
                xP;
                ptsV = vMesh.node(sctrV,1:2);
                pts2 = mesh2.controlPts(sctr,1:2);
                in   = inpolygon(xP(1),xP(2),ptsV(:,1),ptsV(:,2));
                if (in)
                    idu    = mesh2.index(id,1);
                    idv    = mesh2.index(id,2);
                    xiE2   = mesh2.elRangeU(idu,:); % [xi_i,xi_i+1]
                    etaE2  = mesh2.elRangeV(idv,:); % [eta_j,eta_j+1]
                    X2  =  global2LocalMapNURBS2D(xP,pts2,xiE2,etaE2,mesh2.p,mesh2.q,...
                        mesh2.uKnot,mesh2.vKnot,mesh2.weights);
                    GP2 = [GP2;X2 x(3) id];                   
                end

            end
        end
%                         if size(GP2,1) ~= 4 
%                     error('wrong')
%                 end
        couplingSurfaces{is}.quad1{e}=GP1;
        couplingSurfaces{is}.quad2{e}=GP2;
    end
end

%% plot mesh

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
plot_mesh(vMesh.node,vMesh.element(elem2State==1,:),'Q4','g.-',2);
%axis([0 L -c c])

figure
plot_field(vMesh.node,vMesh.element,'Q4',chi);
colorbar

%% find boundary nodes for boundary conditions
% in case of a clamped BC, find next control points

bottomNodes  =  find(abs(mesh2.controlPts(:,2))  <EPS);
topNodes     =  find(abs(mesh2.controlPts(:,2)-L)<EPS);
leftNodes    =  find(abs(mesh2.controlPts(:,1))  <EPS);
rightNodes   =  find(abs(mesh2.controlPts(:,1)-L)<EPS);

fixedNodes  =  unique([bottomNodes;rightNodes;topNodes;leftNodes]);

if clamped
    nextToBotNodes = mesh2.noPtsX+2:2*mesh2.noPtsX-1;
    nextToRgtNodes = 2*mesh2.noPtsX-1:mesh2.noPtsX:mesh2.noPtsX*(mesh2.noPtsY-1)-1;
    nextToTopNodes = mesh2.noPtsX*(mesh2.noPtsY-2)+2:mesh2.noPtsX*(mesh2.noPtsY-1)-1;
    nextToLefNodes = mesh2.noPtsX+2:mesh2.noPtsX:mesh2.noPtsX*(mesh2.noPtsY-2)+2;
    
    nextNodes      = unique([nextToBotNodes';nextToRgtNodes';...
        nextToTopNodes';nextToLefNodes']);
    
    fixedNodes     = [fixedNodes; nextNodes(:)];
end

% constraint inactive plate nodes
fixedNodes = [fixedNodes;inactivePlateNodes];

rightNode = find(node2(:,1)==L);


bndPoints   = node2(rightNode,:);

bndMesh = zeros(mesh2.noElemsV,mesh2.q+1);

for i=1:mesh2.noElemsV
    bndMesh(i,:) = rightNode(i:i+mesh2.q);
end


uFixed    = zeros(1,length(fixedNodes));  % a vector of u_x for the nodes
udofs     = size(mesh1.controlPts,1)*3 + fixedNodes;


numnodes = size(node1,1)   + size(node2,1);
numdofs  = 3*size(node1,1) + size(node2,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE SYSTEM DATA STRUCTURES

disp([num2str(toc),' INITIALIZING DATA STRUCTURES'])

f=zeros(numdofs,1);        % external load vector
K=zeros(numdofs,numdofs);  % stiffness matrix

% ******************************************************************************
% ***                          P R O C E S S I N G                           ***
% ******************************************************************************
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%% COMPUTE STIFFNESS MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp([num2str(toc),'   COMPUTING STIFFNESS MATRIX'])

[W,Q]=quadrature( mesh1.p+1, 'GAUSS', 3 ); % 2x2x2 Gaussian quadrature

%% for domain 1 (continuum)------------------------------------------------

for e=1:size(element1,1)                          % start of element loop
    sctr            = element1(e,:);              % element scatter vector
    nn              = length(sctr);
    sctrB(1:3:3*nn) = 3*sctr-2;
    sctrB(2:3:3*nn) = 3*sctr-1;
    sctrB(3:3:3*nn) = 3*sctr-0;
    sctrf           = 3*sctr-0;
    pts             = node1(sctr,:);
    
    idu    = mesh1.index(e,1);
    idv    = mesh1.index(e,2);
    idw    = mesh1.index(e,3);
    
    xiE    = mesh1.rangeU(idu,:); % [xi_i,xi_i+1]
    etaE   = mesh1.rangeV(idv,:); % [eta_j,eta_j+1]
    zetaE  = mesh1.rangeW(idw,:); % [zeta_k,zeta_k+1]
    
    for q=1:size(W,1)
        pt=Q(q,:);
        wt=W(q);
        
        Xi      = parent2ParametricSpace(xiE,  pt(1));
        Eta     = parent2ParametricSpace(etaE, pt(2));
        Zeta    = parent2ParametricSpace(zetaE,pt(3));
        J2      = jacobianPaPaMapping3d (xiE,etaE,zetaE);
        % compute derivative of basis functions w.r.t parameter coord
        [N dN2dxi dN2deta dRdzeta] = NURBS3DBasisDers([Xi;Eta;Zeta],...
            mesh1.p,mesh1.q,mesh1.r,mesh1.uKnot,mesh1.vKnot,mesh1.wKnot,mesh1.weights');
        jacob      = pts'*[dN2dxi' dN2deta' dRdzeta'];
        J1         = det(jacob);
        dRdx       = [dN2dxi' dN2deta' dRdzeta'] * inv(jacob);
        B          = getBmatrix3D(nn,dRdx);
        K(sctrB,sctrB) = K(sctrB,sctrB) + B' * C * B * J1 * J2 * wt;
        
        f(sctrf)      = f(sctrf)      + q0/t * N' * J1 * J2 * wt;
    end  % of quadrature loop
end

%% for domain 2 (plate)----------------------------------------------------
clear sctrB;

[W,Q]=quadrature(  mesh2.p+1, 'GAUSS', 2 ); % 2 x2point quadrature

for e=1:mesh2.noElems                           % start of element loop
    % skip void elements
    if elem2State(e) == 0
        continue
    end
    sctr   = element2(e,:);                     % element scatter vector
    sctrg  = mesh2.globElems(e,:);              % global element scatter vector
    nn     = length(sctr);
    pts    = node2(sctr,1:2);
    
    idu    = mesh2.index(e,1);
    idv    = mesh2.index(e,2);
    xiE    = mesh2.elRangeU(idu,:); % [xi_i,xi_i+1]
    etaE   = mesh2.elRangeV(idv,:); % [eta_j,eta_j+1]
    
    % loop over Gauss points
    for gp=1:size(W,1)
        pt      = Q(gp,:);
        wt      = W(gp);
        
        % compute coords in parameter space
        Xi      = parent2ParametricSpace(xiE, pt(1));
        Eta     = parent2ParametricSpace(etaE,pt(2));
        J2      = jacobianPaPaMapping(xiE,etaE);
        
        % shape functions, first and second derivatives w.r.t natural coords
        [R dRdxi dRdeta dR2dxi dR2det dR2dxe] = ...
            NURBS2DBasis2ndDers([Xi; Eta],mesh2.p,mesh2.q,mesh2.uKnot,mesh2.vKnot,mesh2.weights');
        
        % compute the jacobian of physical and parameter domain mapping
        % then the derivative w.r.t spatial physical coordinates
        
        jacob  = [dRdxi; dRdeta]          * pts; % 2x2 matrix
        jacob2 = [dR2dxi; dR2det; dR2dxe] * pts; % 3x2 matrix
        
        J1    = det(jacob);
        
        dxdxi = jacob(1,1); dydxi = jacob(1,2);
        dxdet = jacob(2,1); dydet = jacob(2,2);
        
        j33   = [dxdxi^2     dydxi^2     2*dxdxi*dydxi;
            dxdet^2     dydet^2     2*dxdet*dydet;
            dxdxi*dxdet dydxi*dydet dxdxi*dydet+dxdet*dydxi];
        
        dRdx       = inv(jacob)*[dRdxi;dRdeta];
        dR2dx      = inv(j33)  *([dR2dxi; dR2det; dR2dxe]-jacob2*dRdx);
        
        % B matrix
        B              = dR2dx; B(3,:) = B(3,:)*2;
        K(sctrg,sctrg) = K(sctrg,sctrg) + B' * Cp * B * J1 * J2 * wt;
        
        f(sctrg)      = f(sctrg)      + q0 * R' * J1 * J2 * wt;
    end
end

%% interface integrals-----------------------------------------------------

stressState = 'PLANE_STRESS';
% COMPUTE ELASTICITY MATRIX
Cplate = elasticityMatrix(E0,nu0,stressState);
Csolid = C([1 2 4],:);

for is=1:length(couplingSurfaces)
    surf    = couplingSurfaces{is};    
    elemIds = surf.elemIDs;
    quad1   = surf.quad1;
    quad2   = surf.quad2;
                
    for e=1:length(elemIds)
        be       = elemIds (e);        
        sctr1    = element1(be,:);
        pts1     = node1(sctr1,:);        
        nn1      = length(sctr1);                
        sctrB1(1:3:3*nn1) = 3*sctr1-2;
        sctrB1(2:3:3*nn1) = 3*sctr1-1;
        sctrB1(3:3:3*nn1) = 3*sctr1-0;
        
        GP1  = quad1{e};
        GP2  = quad2{e};
                       
        for q=1:size(GP1,1)
            pt1    = GP1(q,1:3);
            pt2    = GP2(q,1:2);
            wt1    = GP1(q,4);
            y      = GP2(q,3);
            normal = GP1(q,5:end);
            n = [normal(1) 0 normal(2);...
                0 normal(2) normal(1);...
                0 0 0];
            n=-n;
            
            e2     = GP2(q,4);
            sctr2  = element2(e2,:);                 % local element scatter vector
            sctr2n = mesh2.globElems(e2,:);
            nn2    = length(sctr2);
            sctrB2 = sctr2n;
            
            pts2 = node2(sctr2,[1 2]);
            
            [N1 dN1dxi dN1deta dN1dzeta] = NURBS3DBasisDers(pt1,...
                mesh1.p,mesh1.q,mesh1.r,mesh1.uKnot,mesh1.vKnot,mesh1.wKnot,mesh1.weights');
            
            J1      = pts1'*[dN1dxi' dN1deta' dN1dzeta'];
            dN1dx   = [dN1dxi' dN1deta' dN1dzeta']*inv(J1);
            B1      = getBmatrix3D(nn1,dN1dx);
            
            [N2 dN2dxi dN2deta dR2dxi dR2det dR2dxe] = ...
                NURBS2DBasis2ndDers(pt2,mesh2.p,mesh2.q,mesh2.uKnot,mesh2.vKnot,mesh2.weights');
            jacob  = [dN2dxi; dN2deta]          * pts2; % 2x2 matrix
            jacob2 = [dR2dxi; dR2det; dR2dxe]   * pts2; % 3x2 matrix
            
            dxdxi = jacob(1,1); dydxi = jacob(1,2);
            dxdet = jacob(2,1); dydet = jacob(2,2);
            
            j33   = [dxdxi^2     dydxi^2     2*dxdxi*dydxi;
                dxdet^2     dydet^2     2*dxdet*dydet;
                dxdxi*dxdet dydxi*dydet dxdxi*dydet+dxdet*dydxi];
            
            % Jacobian inverse and spatial 1st and 2nd derivatives
            
            invJacob   = inv(jacob);
            dRdx       = invJacob*[dN2dxi;dN2deta];
            dN2dx      = inv(j33)*([dR2dxi; dR2det; dR2dxe]-jacob2*dRdx);
            
            % N matrices
            Nm1(1,1:3:3*nn1)  = N1';
            Nm1(2,2:3:3*nn1)  = N1';
            Nm1(3,3:3:3*nn1)  = N1';
            
            Nm2(1,:)  = -y*dRdx(1,:);
            Nm2(2,:)  = -y*dRdx(2,:);
            Nm2(3,:)  =    N2';
            
            B2        = dN2dx;
            B2(3,:)   = B2(3,:)*2;
            B2        = B2*(-y);
            
            % COMPUTE ELEMENT STIFFNESS AT QUADRATURE POINT
            
            Kp11 = alpha*(Nm1'*Nm1)*wt1;
            Kp12 = alpha*(Nm1'*Nm2)*wt1;
            Kp22 = alpha*(Nm2'*Nm2)*wt1;
            
            Kd11 = 0.5 * Nm1'* n * Csolid * B1 *wt1;
            Kd12 = 0.5 * Nm1'* n * Cplate * B2 *wt1;
            Kd21 = 0.5 * Nm2'* n * Csolid * B1 *wt1;
            Kd22 = 0.5 * Nm2'* n * Cplate * B2 *wt1;
            
                    
        K(sctrB1,sctrB1)  = K(sctrB1,sctrB1)  + Kd11 + Kd11' + Kp11;
        K(sctrB1,sctrB2)  = K(sctrB1,sctrB2)  + Kd12 - Kd21' - Kp12;
        K(sctrB2,sctrB1)  = K(sctrB2,sctrB1)  - Kd21 + Kd12' - Kp12';
        K(sctrB2,sctrB2)  = K(sctrB2,sctrB2)  - Kd22 - Kd22' + Kp22;
        
        end  % of quadrature loop
    end
end

K0=K;
%%%%%%%%%%%%%%%%%%% END OF STIFFNESS MATRIX COMPUTATION %%%%%%%%%%%%%%%%%%%


% APPLY ESSENTIAL BOUNDARY CONDITIONS
disp([num2str(toc),'   APPLYING BOUNDARY CONDITIONS'])

bcwt=mean(diag(K)); % a measure of the average  size of an element in K
% used to keep the  conditioning of the K matrix

f=f-K(:,udofs)*uFixed';  % modify the  force vector
f(udofs) = bcwt*uFixed;
K(udofs,:)=0;  % zero out the rows and  columns of the K matrix
K(:,udofs)=0;
K(udofs,udofs)=bcwt*speye(length(udofs));  % put ones*bcwt on the diagonal

% SOLVE SYSTEM
disp([num2str(toc),'   SOLVING SYSTEM'])
U=K\f;

%******************************************************************************
%*** POST - PROCESSING ***
%***************************************************

%%
numnode1 = size(node1,1);
numnode2 = size(node2,1);

index1   = 1:numnode1;
index2   = numnode1+1:numnodes;

Ux1      = U(3*index1-2);
Uy1      = U(3*index1-1);
Uz1      = U(3*index1-0);

Ux2      = U(3*numnode1+1:end);

%% write to VTK
%
% stresses=[sigmaXX sigmaYY sigmaZZ sigmaXY sigmaYZ sigmaZX sigmaVM]

vtuFile1    = 'nitscheIGA-SolidPlate-continuum';
vtuFile2    = 'nitscheIGA-SolidPlate-plate';
vtuFileName = 'nitscheIGA-SolidPlate';

displacement = [Ux1 Uy1 Uz1];
damage       = zeros(length(U),1);

matMap1 = ones(mesh1.noElems,1);
matMap2 = ones(mesh2.noElems,1);
vMesh1=buildVisualizationMesh3D(solid1);
vMesh2=buildVisualizationMesh2D(solid2);
meshes{1}  = mesh1;
meshes{2}  = mesh2;
vmeshes{1} = vMesh1;
vmeshes{2} = vMesh2;
meshData.mesh    = meshes;
meshData.vmesh   = vmeshes;
meshData.matMap{1}=matMap1;
meshData.matMap{2}=matMap2;
material.stiffMat=C;
materials{1}      = material;
meshData.voids   = voidElems;


%figure; hold on;
ok      = plotStress3DForPatch   (meshData,1,vtuFile1,displacement,damage,materials);
ok      = plotStressKirchhoffPlateForPatchEmbed(meshData,2,vtuFile2,Ux2,damage,materials);

pvdFile = fopen(strcat(vtuFileName,'.pvd'), 'wt');

fprintf(pvdFile,'<VTKFile byte_order="LittleEndian" type="Collection" version="0.1">\n');
fprintf(pvdFile,'<Collection>\n');

vtuFile = sprintf('%s%s',vtuFile1,'.vtu');
fprintf(pvdFile,'<DataSet file=''%s'' groups='''' part=''%d'' timestep=''%d''/>\n',...
    vtuFile,0,1);

vtuFile = sprintf('%s%s',vtuFile2,'.vtu');
fprintf(pvdFile,'<DataSet file=''%s'' groups='''' part=''%d'' timestep=''%d''/>\n',...
    vtuFile,0,1);


fprintf(pvdFile,'</Collection>\n');
fprintf(pvdFile,'</VTKFile>\n');

fclose(pvdFile);



