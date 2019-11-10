% This file implements the Nitsche method to join two non-matching meshes.
% Non-matching structured meshes.
%
% Timoshenko beam in bending.
%
% Vinh Phu Nguyen
% Cardiff University, UK
% 1 August 2013

addpath ../../gmshFiles/
addpath ../../post-processing/
addpath ../../fem-functions/
addpath ../../analytical-solutions/
addpath ../../fem_util/;
addpath ../../nurbs-geopdes/inst/
addpath ../../nurbs-util/
addpath ../../meshing/
addpath ../../fem-functions/
addpath ../../post-processing/
addpath ../../xiga/
addpath ../../nurbs-geopdes/inst/
addpath ../../delamination/
addpath ../../xml_toolbox/
addpath ../../finite-cell-method/
addpath ../../C_files/

clear all
clc
state = 0;
tic;

global W Q levelMax

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,'splinecurve.eps',opts)


% MATERIAL PROPERTIES
E0  = 1000;  % Young?s modulus
nu0 = 0.3;   % Poisson?s ratio

E1  = 1;    % Young?s modulus
nu1 = 0.3;  % Poisson?s ratio

% SAMPLE
L  = 5;      % length of the sample
W  = 10;     % width of the sample

gamma = E0/(E1+E0);

plotMesh  = 1;

% FORCE
P = 10; % the peak magnitude of the traction at the right edge

levelMax    = 4;
coord=[-1 -1;1 -1;1 1;-1 1];

% COMPUTE ELASTICITY MATRIX

stressState = 'PLANE_STRESS';
C1 = elasticityMatrix(E0,nu0, stressState);
C2 = elasticityMatrix(E1,nu1, stressState);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE FINITE ELEMENT MESH
%

disp([num2str(toc),'   GENERATING MESH'])

%% mesh1 (background mesh) ------------------------------------------------
controlPts          = zeros(4,2,2);

controlPts(1:2,1,1) = [0;0];
controlPts(1:2,2,1) = [L;0];
controlPts(1:2,1,2) = [0;W];
controlPts(1:2,2,2) = [L;W];

controlPts(4,:,:)   = 1;

uKnot = [0 0 1 1];
vKnot = [0 0 1 1];

solid1 = nrbmak(controlPts,{uKnot vKnot});

% evaluate order

solid1 = nrbdegelev(solid1,[0 0]);

uKnot     = cell2mat(solid1.knots(1));
vKnot     = cell2mat(solid1.knots(2));

refineCountX = 5;
for i=1:refineCountX
    uKnotVectorU = unique(uKnot);
    uKnotVectorV = unique(vKnot);
    newKnotsX = uKnotVectorU(1:end-1) + 0.5*diff(uKnotVectorU);
    newKnotsY = uKnotVectorV(1:end-1) + 0.5*diff(uKnotVectorV);
    newKnots  = {newKnotsX []};
    solid1     = nrbkntins(solid1,newKnots);
    uKnot      = cell2mat(solid1.knots(1));
    vKnot      = cell2mat(solid1.knots(2));
end

refineCountX = 6;
for i=1:refineCountX
    uKnotVectorU = unique(uKnot);
    uKnotVectorV = unique(vKnot);
    newKnotsX = uKnotVectorU(1:end-1) + 0.5*diff(uKnotVectorU);
    newKnotsY = uKnotVectorV(1:end-1) + 0.5*diff(uKnotVectorV);
    newKnots  = {[] newKnotsY};
    solid1     = nrbkntins(solid1,newKnots);
    uKnot      = cell2mat(solid1.knots(1));
    vKnot      = cell2mat(solid1.knots(2));
end

mesh1  = buildIGA2DMesh (solid1);

%% inclusions

r = 1;

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

% move the inclusion to the desired position
controlPts(1,:,:) = controlPts(1,:,:)+L/2;
controlPts(2,:,:) = controlPts(2,:,:)+W/2;

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

% build NURBS object
solid2 = nrbmak(controlPts,{uKnot vKnot});
% p-refinment
solid2 = nrbdegelev(solid2,[0 0]); % to cubic-cubic NURBS
% h-refinement
refineLevel = 4;
for i=1:refineLevel
    uKnotVectorU = unique(uKnot);
    uKnotVectorV = unique(vKnot);
    newKnotsX = uKnotVectorU(1:end-1) + 0.5*diff(uKnotVectorU);
    newKnotsY = uKnotVectorV(1:end-1) + 0.5*diff(uKnotVectorV);
    newKnots  = {newKnotsX newKnotsY};
    solid2     = nrbkntins(solid2,newKnots);
    uKnot      = cell2mat(solid2.knots(1));
    vKnot      = cell2mat(solid2.knots(2));
end

mesh2  = buildIGA2DMesh (solid2);
% do not forget to number mesh2 nodes from mesh1 maximum nodes
mesh2.globElems = mesh2.globElems + size(mesh1.controlPts,1);

%% combine two patches into one data structure

matMap1 = ones(mesh1.elemCount,1);
matMap2 = ones(mesh2.elemCount,1);

vMesh1=buildVisualizationMesh2D(solid1);
vMesh2=buildVisualizationMesh2D(solid2);

meshes{1}  = mesh1;
meshes{2}  = mesh2;
vmeshes{1} = vMesh1;
vmeshes{2} = vMesh2;

meshData.mesh    = meshes;
meshData.vmesh   = vmeshes;
meshData.matMap{1}=matMap1;
meshData.matMap{2}=matMap2;

material1.stiffMat=C1;
material2.stiffMat=C2;

materials{1}      = material1;
materials{2}      = material2;

%% extract boundary curves of solid2 surface
%          (4)
%   -----------------------
%   |                     |
%(1)|                     |(2)
%   |                     |
%   -----------------------
%          (3)

crvs = nrbextract(solid2);

%% boundary mesh for domain 1
bndMesh1 = [];
bndMesh2 = [];


for i=1:mesh2.noElemsV
    bndMesh2 = [bndMesh2; mesh2.noElemsU*i ];
end

for i=1:mesh2.noElemsV
    bndMesh1 = [bndMesh1; mesh2.noElemsU*(i-1)+1 ];
end

bndMesh3 = 1:mesh2.noElemsU;
bndMesh4 = mesh2.noElemsU*(mesh2.noElemsV-1)+1:mesh2.elemCount;

bndMeshes{1}=bndMesh1;
bndMeshes{2}=bndMesh2;
bndMeshes{3}=bndMesh3;
bndMeshes{4}=bndMesh4;

%%

xgp1 = [];
xgp2 = [];
Xgp = [];
wgt = [];
normal=[];


[W1,Q1]=quadrature( 3, 'GAUSS', 1 ); % two point quadrature

for i=1:4
    bnd     = crvs(i);
    bndMesh = bndMeshes{i};
    
    bnd.coefs(1,:) = bnd.coefs(1,:)./bnd.coefs(4,:);
    bnd.coefs(2,:) = bnd.coefs(2,:)./bnd.coefs(4,:);
    
    for e=1:length(bndMesh)
        ptsE    = bnd.coefs(1:2,e:e+mesh2.p)';
        sctr1   = mesh2.locElems(bndMesh(e),:);
        pts1    = mesh2.controlPts(sctr1,:);
        
        if ( i == 1 || i == 2 )
            etE1     = mesh2.rangeV(e,:);   % [et_j,eta_j+1]
        else
            etE1     = mesh2.rangeU(e,:);   % [et_j,eta_j+1]
        end
        
        idu    = mesh2.index(bndMesh(e),1);
        idv    = mesh2.index(bndMesh(e),2);
        xiE    = mesh2.rangeU(idu,:); % [xi_i,xi_i+1]
        etaE   = mesh2.rangeV(idv,:); % [eta_j,eta_j+1]
        
        % order, knot, weights of the NURBS curve bnd
        p       = bnd.order-1;
        knot    = bnd.knots;
        weights = bnd.coefs(4,:);
        
        for q=1:size(W1)
            pt=Q1(q,:);
            wt=W1(q);
            Xi      = 0.5 * ( ( etE1(2) - etE1(1) ) * pt + etE1(2) + etE1(1));
            J2      = 0.5 * ( etE1(2) - etE1(1) );
            [N, dNdxi] = NURBS1DBasisDers(Xi,p,knot,weights);
            J0    = dNdxi*ptsE; % also the tangent to Gamma_*^e
            detJ0 = norm(J0);
            J0    = J0/detJ0;
            
            x     = N*ptsE;
            X1 = global2LocalMapNURBS2D(x,pts1,xiE,etaE,mesh2.p,mesh2.q, ...
                mesh2.uKnot, mesh2.vKnot, mesh2.weights);
            
            Xgp    = [Xgp;x];
            xgp2   = [xgp2;X1 bndMesh(e)];
            wgt    = [wgt;wt*detJ0*J2];
            normal = [normal;-J0(2) J0(1)];            
        end
    end
end

for i=1:length(Xgp)
    x = Xgp(i,:);
   
    
    for id=1:size(vMesh1.element,1)        
        sctrV = vMesh1.element(id,:);
        sctr  = mesh1.locElems(id,:);
        
        ptsV = vMesh1.node(sctrV,1:2);
        pts2 = mesh1.controlPts(sctr,1:2);
        in   = inpolygon(x(1),x(2),ptsV(:,1),ptsV(:,2));
        if (in)
            idu    = mesh1.index(id,1);
            idv    = mesh1.index(id,2);
            xiE2   = mesh1.rangeU(idu,:); % [xi_i,xi_i+1]
            etaE2  = mesh1.rangeV(idv,:); % [eta_j,eta_j+1]
            X2  =  global2LocalMapNURBS2D(x,pts2,xiE2,etaE2,mesh1.p,mesh1.q,...
                                   mesh1.uKnot,mesh1.vKnot,mesh1.weights);
            xgp1 = [xgp1;X2 id];
        end
        
    end
end

%% determination of void, cut and standard elements
% using level set representation of  the coupling surface


nNode = size(vMesh1.node,1);
chi   = zeros(1,nNode);

xc = L/2; % x coord of center
yc = W/2; % y coord of center

for i = 1 : nNode
    x      = vMesh1.node(i,1);
    y      = vMesh1.node(i,2);
    d      = sqrt((x-xc)^2+(y-yc)^2);
    chi(i) = d - r; % level set
end

elem1State = zeros(size(mesh1.locElems,1),1);

% loop over elements

for iel = 1 : size(mesh1.locElems,1)    
    sctrV   = vMesh1.element(iel,:);
    
    phi   = chi(1,sctrV);
    maPhi = max(phi);
    miPhi = min(phi);
    
    if     ( maPhi*miPhi < 0 )
        abs(maPhi/miPhi)
        elem1State(iel) = 1;             % cut elements
    elseif ( maPhi <= 0 )
        elem1State(iel) = 0;             % void elements
    else
        elem1State(iel) = 2;             % normal elements
    end
    
    if abs(maPhi/miPhi) > 10 
        elem1State(iel) = 2;             % void elements
    end
    
    if abs(maPhi/miPhi) < 0.003
        elem1State(iel) = 0;             % void elements
    end
end

% find inactive plate dofs

inactiveNodes1 = [];
inactiveNodes2 = [];

for i=1:size(mesh1.locElems,1)
    sctr = mesh1.locElems(i,:);
    if     elem1State(i) == 0
        inactiveNodes1 = [inactiveNodes1;sctr];
    elseif elem1State(i) == 1
        inactiveNodes2 = [inactiveNodes2;sctr];
    end
end

inactiveNodes1 = unique(inactiveNodes1);
inactiveNodes2 = unique(inactiveNodes2);

inactiveNodes = setdiff(inactiveNodes1,inactiveNodes2);

voidElems = find(elem1State==0);
cutElems  = find(elem1State==1);

meshData.voids   = voidElems;

%% DEFINE BOUNDARIES

% GET NODES ON DISPLACEMENT BOUNDARY
%      Here we get the nodes on the essential boundaries

fixedNode=find(mesh1.controlPts(:,2)==0);
topNode  =find(mesh1.controlPts(:,2)==W);

topEdge   = zeros(mesh1.noElemsU,mesh1.p+1);

for i=1:mesh1.noElemsU
    topEdge(i,:) = topNode(i:i+mesh1.p);
end

fixedNode = [fixedNode;inactiveNodes];

uFixed=zeros(1,length(fixedNode))';  % a vector of u_x for the nodes
vFixed=zeros(1,length(fixedNode))';



%% plot meshes

figure
hold on
nrbkntplot(solid1);
nrbkntplot(solid2);
 %plot(mesh1.controlPts(:,1),mesh1.controlPts(:,2),'o','MarkerEdgeColor','k',...
 %    'MarkerFaceColor','g','MarkerSize',6.5);
 %plot(mesh2.controlPts(:,1),mesh2.controlPts(:,2),'s','MarkerEdgeColor','k',...
 %    'MarkerFaceColor','r','MarkerSize',6.5);
plot(Xgp(:,1),Xgp(:,2),'s','MarkerEdgeColor','k',...
    'MarkerFaceColor','r','MarkerSize',6.5);
plot_mesh(vMesh1.node,vMesh1.element(cutElems,:),'Q4','r-',1.2);
plot_mesh(vMesh1.node,vMesh1.element(voidElems,:),'Q4','cy-',1.4);
axis off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE SYSTEM DATA STRUCTURES

disp([num2str(toc),' INITIALIZING DATA STRUCTURES'])

numnodes = size(mesh1.controlPts,1) + size(mesh2.controlPts,1);

f=zeros(2*numnodes,1);          % external load vector
K=zeros(2*numnodes,2*numnodes); % stiffness matrix

% ******************************************************************************
% ***                          P R O C E S S I N G                           ***
% ******************************************************************************
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%% COMPUTE STIFFNESS MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp([num2str(toc),'   COMPUTING STIFFNESS MATRIX'])


%% for domain 1 -----------------------------------------------------------

for ip=1:length(meshes)
    mesh = meshes{ip};
    [W,Q]= quadrature( mesh.p+1, 'GAUSS', 2 ); % 2x2 Gaussian quadrature
    C = materials{ip}.stiffMat;
    
    for e=1:mesh.elemCount
        % if background patch and void element then skip it
        if (ip ==1 && elem1State(e) == 0)
            continue
        end
        % if background patch and cut element then special Gauss quad
        if (ip ==1 && ismember(e,cutElems))     
            W = [];
            Q = [];
            [aa] = hierarchicalGaussQuad(mesh.p+1,chi(vMesh1.element(e,:)),coord,0);
            CN        = chi(vMesh1.element(e,:));
            xi        = Q(:,1);
            eta       = Q(:,2);
            N         = 1/4*[(1-xi).*(1-eta) (1+xi).*(1-eta) (1+xi).*(1+eta) (1-xi).*(1+eta)];
            ls        = N(:,1)*CN(1)+N(:,2)*CN(2)+N(:,3)*CN(3)+N(:,4)*CN(4);
            id         = find(ls <= 0);
            W(id)      = [];
            Q(id,:)    = [];
        else
            [W,Q]= quadrature( mesh.p+1, 'GAUSS', 2 ); 
        end
        
        idu    = mesh.index(e,1);
        idv    = mesh.index(e,2);
        xiE    = mesh.rangeU(idu,:); % [xi_i,xi_i+1]
        etaE   = mesh.rangeV(idv,:); % [eta_j,eta_j+1]
        
        sctrg  = mesh.globElems(e,:);         %  element scatter vector
        sctrl  = mesh.locElems (e,:);         %  element scatter vector
        nn     = length(sctrg);
        %sctrB  = zeros(1,2*nn);
        sctrB(1,1:2:2*nn) = 2*sctrg-1;
        sctrB(1,2:2:2*nn) = 2*sctrg  ;
        
        pts    = mesh.controlPts(sctrl,:);
        
        % loop over Gauss points
        for gp=1:size(W,1)
            pt      = Q(gp,:);
            wt      = W(gp);
            Xi      = parent2ParametricSpace(xiE,pt(1));
            Eta     = parent2ParametricSpace(etaE,pt(2));
            J2      = jacobianPaPaMapping(xiE,etaE);
            
            % compute derivatives of shape functions
            [R, dRdxi, dRdeta] = NURBS2DBasisDers([Xi;Eta],mesh.p,mesh.q,mesh.uKnot,mesh.vKnot,mesh.weights');
            
            % compute the jacobian of physical and parameter domain mapping
            % then the derivative w.r.t spatial physical coordinates
            
            jacob      =  pts' * [dRdxi' dRdeta'];
            J1         = det(jacob);
            invJacob   = inv(jacob);
            dRdx       = [dRdxi' dRdeta'] * invJacob;
            J          = J1 * J2;
            B           = getBmatrix2D(dRdx');
            K(sctrB,sctrB) = K(sctrB,sctrB) + B' * C * B * J * wt;
        end
    end
end


%% interface integrals ----------------------------------------------------

% stabilisation param

beta = 1e6;

for i=1:length(wgt)                     % start of element loop
    e1     = xgp1(i,3);
    e2     = xgp2(i,3);
    sctr1  = mesh1.globElems(e1,:);
    sctr2  = mesh2.globElems(e2,:);
    sctr2L = mesh2.locElems (e2,:);
    
    nn1 = length(sctr1);
    nn2 = length(sctr2);
    
    sctrB1(1,1:2:2*nn1) = 2*sctr1-1;
    sctrB1(1,2:2:2*nn1) = 2*sctr1  ;
    
    sctrB2(1,1:2:2*nn2) = 2*sctr2-1;
    sctrB2(1,2:2:2*nn2) = 2*sctr2  ;
    
    pts1 = mesh1.controlPts(sctr1,:);
    pts2 = mesh2.controlPts(sctr2L,:);
    
    idv1   = mesh1.index(e1,2);
    idv2   = mesh2.index(e2,2);
    xiE1      = mesh1.rangeU(end,:); % [xi_i,xi_i+1]
    etE1      = mesh1.rangeV(idv1,:);   % [et_j,eta_j+1]
    xiE2      = mesh2.rangeU(1,:); % [xi_i,xi_i+1]
    etE2      = mesh2.rangeV(idv2,:);   % [et_j,eta_j+1]
    
    pt1 = xgp1(i,1:2);
    pt2 = xgp2(i,1:2);
    wt1 = wgt(i);
    
    normalv=normal(i,:);
    n = [normalv(1) 0 normalv(2);0 normalv(2) normalv(1)];
    n = -n;
    
    % compute derivatives of shape functions
    [R1 dR1dxi dR1deta] = NURBS2DBasisDers(pt1,mesh1.p,mesh1.q,mesh1.uKnot,mesh1.vKnot,mesh1.weights');
    [R2 dR2dxi dR2deta] = NURBS2DBasisDers(pt2,mesh2.p,mesh2.q,mesh2.uKnot,mesh2.vKnot,mesh2.weights');
    
    % compute the jacobian of physical and parameter domain mapping
    % then the derivative w.r.t spatial physical coordinates
    
    jacob1     = [dR1dxi; dR1deta] * pts1;
    jacob2     = [dR2dxi; dR2deta] * pts2;
    invJacob1   = inv(jacob1);
    invJacob2   = inv(jacob2);
    dR1dx      = [dR1dxi' dR1deta'] * invJacob1;
    dR2dx      = [dR2dxi' dR2deta'] * invJacob2;
    
    
    B1(1,1:2:nn1*2)  = dR1dx(:,1)';
    B1(2,2:2:nn1*2)  = dR1dx(:,2)';
    B1(3,1:2:nn1*2)  = dR1dx(:,2)';
    B1(3,2:2:nn1*2)  = dR1dx(:,1)';
    
    B2(1,1:2:nn2*2)  = dR2dx(:,1)';
    B2(2,2:2:nn2*2)  = dR2dx(:,2)';
    B2(3,1:2:nn2*2)  = dR2dx(:,2)';
    B2(3,2:2:nn2*2)  = dR2dx(:,1)';
    
    Nm1(1,1:2:nn1*2)  = R1';
    Nm1(2,2:2:nn1*2)  = R1';
    Nm2(1,1:2:nn2*2)  = R2';
    Nm2(2,2:2:nn2*2)  = R2';
    
    % COMPUTE ELEMENT STIFFNESS AT QUADRATURE POINT
    
    Kp11 = beta*(Nm1'*Nm1)*wt1;
    Kp12 = beta*(Nm1'*Nm2)*wt1;
    Kp22 = beta*(Nm2'*Nm2)*wt1;
    
    Kd11 = gamma     * Nm1'* n * C1 * B1 *wt1;
    Kd12 = (1-gamma) * Nm1'* n * C2 * B2 *wt1;
    Kd21 = gamma     * Nm2'* n * C1 * B1 *wt1;
    Kd22 = (1-gamma) * Nm2'* n * C2 * B2 *wt1;
    
    
    K(sctrB1,sctrB1)  = K(sctrB1,sctrB1)  - Kd11 - Kd11' + Kp11;
    K(sctrB1,sctrB2)  = K(sctrB1,sctrB2)  - Kd12 + Kd21' - Kp12;
    K(sctrB2,sctrB1)  = K(sctrB2,sctrB1)  + Kd21 - Kd12' - Kp12';
    K(sctrB2,sctrB2)  = K(sctrB2,sctrB2)  + Kd22 + Kd22' + Kp22;
end



%% External force

[W,Q]=quadrature( 3, 'GAUSS', 1 ); % three point quadrature

for e=1:mesh1.noElemsU
    xiE   = mesh1.rangeU(e,:); % [xi_i,xi_i+1]
    conn  = topEdge(e,:);
    pts   = mesh1.controlPts(conn,:);
    sctr  = topEdge(e,:);
    sctrx = 2*sctr-1;
    sctry = 2*sctr;
    
    % loop over Gauss points
    for gp=1:size(W1,1)
        xi      = Q1(gp,:);
        wt      = W1(gp);
        Xi      = 0.5 * ( ( xiE(2) - xiE(1) ) * xi + xiE(2) + xiE(1));
        J2      = 0.5 * ( xiE(2) - xiE(1) );
        
        [N, dNdxi] = NURBS1DBasisDers(Xi,mesh1.p,mesh1.uKnot,mesh1.weights);
        
        % compute the jacobian of physical and parameter domain mapping
        % then the derivative w.r.t spatial physical coordinates
        
        x        = N     * pts; % global coord of GP
        jacob1   = dNdxi * pts;
        J1       = norm (jacob1);
        fac = J1 * J2 * wt;  
        f(sctry) = f(sctry) + N' * P * fac;
    end
end


%%%%%%%%%%%%%%%%%%% END OF STIFFNESS MATRIX COMPUTATION %%%%%%%%%%%%%%%%%%%


%% APPLY ESSENTIAL BOUNDARY CONDITIONS
disp([num2str(toc),'   APPLYING BOUNDARY CONDITIONS'])
bcwt=mean(diag(K)); % a measure of the average size of an element in K
% used to keep the conditioning of the K matrix
udofs=2*fixedNode-1;           % global indecies of the fixed x displacements
vdofs=2*fixedNode;   % global indecies of the fixed y displacements
f=f-K(:,udofs)*uFixed;  % modify the force vector
f=f-K(:,vdofs)*vFixed;
K(udofs,:)=0;
K(vdofs,:)=0;
K(:,udofs)=0;
K(:,vdofs)=0;
K(udofs,udofs)=bcwt*speye(length(udofs)); % put ones*bcwt on the diagonal
K(vdofs,vdofs)=bcwt*speye(length(vdofs));
f(udofs)=bcwt*speye(length(udofs))*uFixed;
f(vdofs)=bcwt*speye(length(udofs))*vFixed;

% SOLVE SYSTEM
disp([num2str(toc),'   SOLVING SYSTEM'])
U=K\f;

%**************************************************************************
%*** POST - PROCESSING ***
%***************************************************

disp([num2str(toc),'   POST-PROCESSING'])

%% write to Paraview
displacement = [U(1:2:end) U(2:2:end)];


vtuFileName = 'iga-nitsche-inclusion';
for ip=1:2
    vtuFile = strcat(vtuFileName,'-mesh',num2str(ip));
    %figure; hold on;
    ok      = plotStress2DForPatch(meshData,ip,vtuFile,displacement,materials);
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


