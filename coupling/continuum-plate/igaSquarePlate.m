%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Isogeometric analysis for Kirchoff plate problems.
%
% Rotation-free thin plates. Fully clamped rectangular plates.
% To compare solid-plate coupled model with.
%
% Vinh Phu Nguyen,
% Cardiff University
% nvinhphu@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath ../fem_util/
addpath ../C_files/
addpath ../data/
addpath ../meshing/
addpath ../post-processing/
addpath ../fem-functions/
addpath ../meshing/
addpath ../nurbs-util/
addpath ../nurbs-geopdes/inst/

clc
clear all

tic;


% MATERIAL PROPERTIES
E  = 1000;  % Young modulus
nu = 0.3;   % Poisson ratio

% Plate PROPERTIES
L     = 400;     % length of the beam
t     = 20;      % thicknes

q0    = -10;    % uniformly distributed loads


% constitutive matrix

D  = E*t^3/(12*(1-nu^2));
C  = D*[1 nu 0;nu 1 0; 0 0 0.5*(1-nu)];

%% mesh

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

mesh     = buildIGA2DMeshForSurface (solid2);

%% find boundary nodes for boundary conditions

EPS = 1e-8;
bottomNodes  =  find(abs(mesh.controlPts(:,2))  <EPS);
topNodes     =  find(abs(mesh.controlPts(:,2)-L)<EPS);
leftNodes    =  find(abs(mesh.controlPts(:,1))  <EPS);
rightNodes   =  find(abs(mesh.controlPts(:,1)-L)<EPS);

fixedNodes  =  unique([bottomNodes;rightNodes;topNodes;leftNodes]);

clamped=1;

if clamped
nextToBotNodes = mesh.noPtsX+2:2*mesh.noPtsX-1;
nextToRgtNodes = 2*mesh.noPtsX-1:mesh.noPtsX:mesh.noPtsX*(mesh.noPtsY-1)-1;
nextToTopNodes = mesh.noPtsX*(mesh.noPtsY-2)+2:mesh.noPtsX*(mesh.noPtsY-1)-1;
nextToLefNodes = mesh.noPtsX+2:mesh.noPtsX:mesh.noPtsX*(mesh.noPtsY-2)+2;

nextNodes      = unique([nextToBotNodes';nextToRgtNodes';...
                         nextToTopNodes';nextToLefNodes']);

fixedNodes     = [fixedNodes; nextNodes(:)];
end

plot(mesh.controlPts(fixedNodes,1),mesh.controlPts(fixedNodes,2),...
    'bs','MarkerEdgeColor','r','MarkerSize',14);


noCtrPts       = mesh.noPts;
noDofs         = noCtrPts * 1;

%% initialization

K = sparse(noDofs,noDofs);  % global stiffness matrix
f = zeros(noDofs,1);        % external force vector

% essential boundary conditions

uFixed     = zeros(size(fixedNodes))';
udofs      = fixedNodes;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gauss quadrature rule
noGPs = mesh.p+1;
noGpEle = noGPs^2;
[W,Q]=quadrature(  noGPs, 'GAUSS', 2 ); % noGPs x noGPs point quadrature

% Assembling system of equation
% Stiffness matrix and external force vector

disp([num2str(toc),'  ASSEMBLING THE SYSTEM']);

% Loop over elements (knot spans)

for e=1:mesh.noElems
    idu    = mesh.index(e,1);
    idv    = mesh.index(e,2);
    xiE    = mesh.elRangeU(idu,:); % [xi_i,xi_i+1]
    etaE   = mesh.elRangeV(idv,:); % [eta_j,eta_j+1]
    
    sctr   = mesh.locElems(e,:);          %  element scatter vector
    nn     = length(sctr);
    pts    = mesh.controlPts(sctr,1:2);
    
    % loop over Gauss points
    
    for gp=1:size(W,1)
        pt      = Q(gp,:);
        wt      = W(gp);
        
        % compute coords in parameter space
        Xi      = parent2ParametricSpace(xiE,pt(1));
        Eta     = parent2ParametricSpace(etaE,pt(2));
        J2      = jacobianPaPaMapping(xiE,etaE);
        
        % shape functions, first and second derivatives w.r.t natural coords
        
        [R dRdxi dRdeta dR2dxi dR2det dR2dxe] = ...
            NURBS2DBasis2ndDers([Xi;Eta],mesh.p,mesh.q,mesh.uKnot,mesh.vKnot,mesh.weights');
        
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
        
        % Jacobian inverse and spatial 1st and 2nd derivatives
        
        invJacob   = inv(jacob);
        dRdx       = invJacob*[dRdxi;dRdeta];
        dR2dx      = inv(j33)*([dR2dxi; dR2det; dR2dxe]-jacob2*dRdx);
        
        % B matrix
        
        B          = dR2dx;
        B(3,:)     = B(3,:)*2;
        
        % compute elementary stiffness matrix and
        % assemble it to the global matrix
        
        K(sctr,sctr) = K(sctr,sctr) + B' * C * B * J1 * J2 * wt;
        f(sctr)      = f(sctr)      + q0 * R' * J1 * J2 * wt;
    end
end

disp([num2str(toc),'  APPLYING BOUNDARY  CONDITIONS']);
bcwt=mean(diag(K)); % a measure of the average  size of an element in K
% used to keep the  conditioning of the K matrix

f=f-K(:,udofs)*uFixed';  % modify the  force vector
f(udofs) = bcwt*uFixed;
K(udofs,:)=0;  % zero out the rows and  columns of the K matrix
K(:,udofs)=0;
K(udofs,udofs)=bcwt*speye(length(udofs));  % put ones*bcwt on the diagonal

% SOLVE SYSTEM
disp([num2str(toc),'  SOLVING THE SYSTEM']);
U=K\f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% POST-PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Post-processing using triangulation of GPs

%% write to VTK
%
% stresses=[sigmaXX sigmaYY sigmaZZ sigmaXY sigmaYZ sigmaZX sigmaVM]

vtuFile    = 'squarePlate';



damage       = zeros(length(U),1);

matMap1 = ones(mesh.noElems,1);
vMesh=buildVisualizationMesh2D(solid2);
meshes{1}  = mesh;
vmeshes{1} = vMesh;
meshData.mesh    = meshes;
meshData.vmesh   = vmeshes;
meshData.matMap{1}=matMap1;
material.stiffMat=C;
materials{1}      = material;

%figure; hold on;

ok      = plotStressKirchhoffPlateForPatch(meshData,1,vtuFile,U,damage,materials);






