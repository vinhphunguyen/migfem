% This file implements the Nitsche method to join two non-matching meshes.
% Non-matching meshes version. NURBS elements.
% Non hierarchical meshes and one big patch is coupled to a smaller patch.
%
% Timoshenko beam.
%
% -----------------------
% |                     |-------------------
% |                     |                  |
% |                     |                  |
% |                     |-------------------
% -----------------------
%
%
% Vinh Phu Nguyen
% Cardiff University, UK
% 26 July 2013

addpath ../../gmshFiles/
addpath ../../post-processing/
addpath ../../fem-functions/
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
addpath ../../C_files/

clear all
colordef white
state = 0;
tic;

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,'splinecurve.eps',opts)

% penalty (reguralisation) parameter

alpha = 1e9;

symmetric = 1;

if (symmetric)
    %Symmetric Nitsche method
    gamma1 = 1;
    gamma2 = 1;
else
    %Non-Symmetric Nitsche method
    gamma1 = -1;
    gamma2 = 0;
end

% MATERIAL PROPERTIES
E0  = 3e7;  % Young?s modulus
nu0 = 0.3;  % Poisson?s ratio

% BEAM PROPERTIES
L  = 48;     % length of the beam
W1 = 6;
W2 = 4;
c1 = W1/2;
c2 = W2/2;

L1=L/2;


% TIP LOAD
F = -1000; % the peak magnitude of the traction at the right edge

% COMPUTE ELASTICITY MATRIX
C=E0/(1-nu0^2)*[   1      nu0          0;
    nu0        1          0;
    0        0  (1-nu0)/2 ];

material.stiffMat=C;

%% MESHING

%% mesh1 ------------------------------------------------------------------

controlPts          = zeros(4,2,2);

controlPts(1:2,1,1) = [0;-c1];
controlPts(1:2,2,1) = [L1;-c1];
controlPts(1:2,1,2) = [0;c1];
controlPts(1:2,2,2) = [L1;c1];

controlPts(4,:,:)   = 1;

uKnot = [0 0 1 1];
vKnot = [0 0 1 1];

solid1 = nrbmak(controlPts,{uKnot vKnot});

% evaluate order

solid1 = nrbdegelev(solid1,[2 2]);

uKnot     = cell2mat(solid1.knots(1));
vKnot     = cell2mat(solid1.knots(2));

%solid1     = nrbkntins(solid1,{[] [0.3 0.7]});

refineCountX = 6;
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

refineCountX = 4;
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

mesh1 = buildIGA2DMesh (solid1);

%% mesh2 ------------------------------------------------------------------

controlPts          = zeros(4,2,2);

controlPts(1:2,1,1) = [L1;-c2];
controlPts(1:2,2,1) = [L;-c2];
controlPts(1:2,1,2) = [L1;c2];
controlPts(1:2,2,2) = [L;c2];

controlPts(4,:,:)   = 1;

uKnot = [0 0 1 1];
vKnot = [0 0 1 1];

solid2 = nrbmak(controlPts,{uKnot vKnot});

% evaluate order

solid2 = nrbdegelev(solid2,[2 2]);

uKnot     = cell2mat(solid2.knots(1));
vKnot     = cell2mat(solid2.knots(2));

refineCountX = 6;
for i=1:refineCountX
    uKnotVectorU = unique(uKnot);
    uKnotVectorV = unique(vKnot);
    newKnotsY = [];
    newKnotsX = uKnotVectorU(1:end-1) + 0.5*diff(uKnotVectorU);
    newKnotsY = uKnotVectorV(1:end-1) + 0.5*diff(uKnotVectorV);
    newKnots  = {newKnotsX []};
    solid2     = nrbkntins(solid2,newKnots);
    uKnot      = cell2mat(solid2.knots(1));
    vKnot      = cell2mat(solid2.knots(2));
end

%solid2     = nrbkntins(solid2,{[] [0.2 0.4 0.6 0.8]});

refineCountX = 2;
for i=1:refineCountX
    uKnotVectorU = unique(uKnot);
    uKnotVectorV = unique(vKnot);
    newKnotsY = [];
    newKnotsX = uKnotVectorU(1:end-1) + 0.5*diff(uKnotVectorU);
    newKnotsY = uKnotVectorV(1:end-1) + 0.5*diff(uKnotVectorV);
    newKnots  = {[] newKnotsY};
    solid2     = nrbkntins(solid2,newKnots);
    uKnot      = cell2mat(solid2.knots(1));
    vKnot      = cell2mat(solid2.knots(2));
end

mesh2 = buildIGA2DMesh (solid2);
mesh2.globElems = mesh2.globElems + size(mesh1.controlPts,1);

% plot meshes

figure
hold on
nrbkntplot(solid1);
nrbkntplot(solid2);
plot(mesh1.controlPts(:,1),mesh1.controlPts(:,2),'o','MarkerEdgeColor','k',...
    'MarkerFaceColor','g','MarkerSize',4.5);
plot(mesh2.controlPts(:,1),mesh2.controlPts(:,2),'s','MarkerEdgeColor','k',...
    'MarkerFaceColor','g','MarkerSize',6.5);
axis off

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
materials{1}      = material;

%% boundary mesh for domain 1
bndMesh1 = [];
for i=1:mesh1.noElemsV
    bndMesh1 = [bndMesh1; mesh1.noElemsU*i ];
end

%% boundary mesh for domain 2
bndMesh2 = [];
for i=1:mesh2.noElemsV
    bndMesh2 = [bndMesh2; mesh2.noElemsU*(i-1)+1 ];
end

%% boundary edges (to distribute GPs)

bndNodes2  = find(abs(mesh2.controlPts(:,1)-L1)<1e-14);
bndEdge2   = zeros(mesh2.noElemsV,mesh2.q+1);

for i=1:mesh2.noElemsV
    bndEdge2(i,:) = bndNodes2(i:i+mesh2.q);
end

bndNodes1  = find(abs(vMesh1.node(:,1)-L1)<1e-14);
bndCoords1 = vMesh1.node(bndNodes1,2);

%% determine GPs (for mesh1 and mesh2) on coupling interface

ngp = (mesh1.p+1)*(mesh1.q+1);

[W1,Q1]=quadrature( mesh2.q+1, 'GAUSS', 1 ); % two point quadrature

xgp1 = [];
xgp2 = [];
Xgp = [];
wgt = [];
normal=[];
elems1=[];
elems2=[];

for e=1:mesh2.noElemsV
    sctrEdge = bndEdge2(e,:);
    sctr2    = mesh2.locElems(bndMesh2(e),:);
    pts2     = mesh2.controlPts(sctr2,:);
    ptsS2    = mesh2.controlPts(sctrEdge,:);
    xiE2     = mesh2.rangeU(1,:);   % [xi_i,xi_i+1]
    etE2     = mesh2.rangeV(e,:);   % [et_j,eta_j+1]
    for q=1:size(W1)
        pt=Q1(q,:);
        wt=W1(q);
        Xi      = 0.5 * ( ( etE2(2) - etE2(1) ) * pt + etE2(2) + etE2(1));
        J2      = 0.5 *   ( etE2(2) - etE2(1) );
        [N dNdxi] = NURBS1DBasisDers(Xi,mesh2.q,mesh2.vKnot,mesh2.weights);
        J0        = dNdxi*ptsS2;          % also the tangent to Gamma_*^e
        detJ0     = norm(J0);
        J0        = J0/detJ0;
        
        x         = N*ptsS2; % global coords of GP q
        aa        = bndCoords1-x(2);
        bb        = find(aa>0);
        id        = bb(1) - 1;
        
        X2 = global2LocalMapNURBS2D(x,pts2,xiE2,etE2,mesh2.p,mesh2.q, ...
            mesh2.uKnot, mesh2.vKnot, mesh2.weights);
        
        Xgp    = [Xgp;x];
        xgp2   = [xgp2;X2];
        wgt    = [wgt;wt*detJ0*J2];
        normal = [normal;-J0(2) J0(1)];        
        elems2 = [elems2;bndMesh2(e)];        
        elems1 = [elems1;bndMesh1(id)];
        
        e1     = bndMesh1(id);
        idu    = mesh1.index(e1,1);
        idv    = mesh1.index(e1,2);
        sctr1    = mesh1.locElems(e1,:);
        pts1     = mesh1.controlPts(sctr1,:);
        xiE1     = mesh1.rangeU(idu,:); % [xi_i,xi_i+1]
        etE1     = mesh1.rangeV(idv,:);   % [et_j,eta_j+1]
        X1 = global2LocalMapNURBS2D(x,pts1,xiE1,etE1,mesh1.p,mesh1.q, ...
            mesh1.uKnot, mesh1.vKnot, mesh1.weights);
        xgp1 = [xgp1;X1];
    
    end
end


% DEFINE BOUNDARIES

% GET NODES ON DISPLACEMENT BOUNDARY
%      Here we get the nodes on the essential boundaries

fixedNode=find(mesh1.controlPts(:,1)==0);
rightNode=find(mesh2.controlPts(:,1)==L);

rightEdge   = zeros(mesh2.noElemsV,mesh2.q+1);

for i=1:mesh2.noElemsV
    rightEdge(i,:) = rightNode(i:i+mesh2.q);
end

uFixed=zeros(1,length(fixedNode))';  % a vector of u_x for the nodes
vFixed=zeros(1,length(fixedNode))';

numnodes = size(mesh1.controlPts,1) + size(mesh2.controlPts,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE SYSTEM DATA STRUCTURES

disp([num2str(toc),' INITIALIZING DATA STRUCTURES'])

f=zeros(2*numnodes,1);          % external load vector
K=zeros(2*numnodes,2*numnodes); % stiffness matrix

% ******************************************************************************
% ***                          P R O C E S S I N G                           ***
% ******************************************************************************
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%% COMPUTE STIFFNESS MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp([num2str(toc),'   COMPUTING STIFFNESS MATRIX'])

%% continuum part----------------------------------------------------------

for ip=1:length(meshes)
    mesh = meshes{ip};
    [W,Q]= quadrature( mesh.p+1, 'GAUSS', 2 ); % 2x2 Gaussian quadrature
    
    for e=1:mesh.elemCount
        idu    = mesh.index(e,1);
        idv    = mesh.index(e,2);
        xiE    = mesh.rangeU(idu,:); % [xi_i,xi_i+1]
        etaE   = mesh.rangeV(idv,:); % [eta_j,eta_j+1]
        
        sctrg  = mesh.globElems(e,:);         %  element scatter vector
        sctrl  = mesh.locElems (e,:);         %  element scatter vector
        nn     = length(sctrg);
        sctrB  = zeros(1,2*nn);
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
            [R dRdxi dRdeta] = NURBS2DBasisDers([Xi; Eta],mesh.p,mesh.q,mesh.uKnot,mesh.vKnot,mesh.weights');
            
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

for i=1:length(wgt)                     % start of element loop
    e1     = elems1(i);
    e2     = elems2(i);
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
    
    Kp11 = alpha*(Nm1'*Nm1)*wt1;
    Kp12 = alpha*(Nm1'*Nm2)*wt1;
    Kp22 = alpha*(Nm2'*Nm2)*wt1;
    
    Kd11 = 0.5 * Nm1'* n * C * B1 *wt1;
    Kd12 = 0.5 * Nm1'* n * C * B2 *wt1;
    Kd21 = 0.5 * Nm2'* n * C * B1 *wt1;
    Kd22 = 0.5 * Nm2'* n * C * B2 *wt1;
    

    K(sctrB1,sctrB1)  = K(sctrB1,sctrB1)  - Kd11 - gamma1*Kd11' + gamma2*Kp11;
    K(sctrB1,sctrB2)  = K(sctrB1,sctrB2)  - Kd12 + gamma1*Kd21' - gamma2*Kp12;
    K(sctrB2,sctrB1)  = K(sctrB2,sctrB1)  + Kd21 - gamma1*Kd12' - gamma2*Kp12';
    K(sctrB2,sctrB2)  = K(sctrB2,sctrB2)  + Kd22 + gamma1*Kd22' + gamma2*Kp22;
end

%% External force ---------------------------------------------------------

[W,Q]=quadrature( 9, 'GAUSS', 1 ); % three point quadrature

for e=1:mesh2.noElemsV
    xiE   = mesh2.rangeV(e,:); % [xi_i,xi_i+1]
    conn  = rightEdge(e,:);
    pts   = mesh2.controlPts(conn,:);
    sctr  = rightEdge(e,:) + size(mesh1.controlPts,1);
    sctrx = 2*sctr-1;
    sctry = 2*sctr;
    
    % loop over Gauss points
    for gp=1:size(W1,1)
        xi      = Q1(gp,:);
        wt      = W1(gp);
        Xi      = 0.5 * ( ( xiE(2) - xiE(1) ) * xi + xiE(2) + xiE(1));
        J2      = 0.5 * ( xiE(2) - xiE(1) );
        
        [N dNdxi] = NURBS1DBasisDers(Xi,mesh2.q,mesh2.vKnot,mesh2.weights);
        
        % compute the jacobian of physical and parameter domain mapping
        % then the derivative w.r.t spatial physical coordinates
        
        x        = N     * pts; % global coord of GP
        jacob1   = dNdxi * pts;
        J1       = norm (jacob1);
        fac      = J1 * J2 * wt;
        f(sctry) = f(sctry) + N' * F * fac;
    end
end


%%%%%%%%%%%%%%%%%%% END OF STIFFNESS MATRIX COMPUTATION %%%%%%%%%%%%%%%%%%%


% APPLY ESSENTIAL BOUNDARY CONDITIONS
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

%******************************************************************************
%*** POST - PROCESSING ***
%***************************************************

disp([num2str(toc),'   POST-PROCESSING'])

displacement = [U(1:2:end) U(2:2:end)];


vtuFileName = 'iga-nitsche-tbeam-general-dim';
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
%%

