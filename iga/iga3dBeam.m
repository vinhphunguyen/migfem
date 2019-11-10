%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Isogeometric analysis for three dimensional elasticity problems.
%
% 3D beam in bending
%
% Vinh Phu Nguyen,
% Johns Hopkins University
% nvinhphu@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath ../fem_util/
addpath ../C_files/
addpath ../data/
addpath ../meshing/
addpath ../integration/
addpath ../post-processing/
addpath ../fem-functions/
addpath ../nurbs-util/
addpath ../analytical-solutions/

clc
clear all

global p q r

E0           = 1e5;  % Young modulus
nu0          = 0.3;  % Poisson’s ratio

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE COMPLIANCE MATRIX
D=zeros(6,6);
D(1:3,1:3)=E0/(1+nu0)/(1-2*nu0)*[ 1-nu0 nu0 nu0;
                                  nu0 1-nu0 nu0;
                                  nu0 nu0 1-nu0];
D(4:6,4:6)=E0/2/(1+nu0)*eye(3);

tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%  CAD input: control points and knot vectors
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

beam3dC1Data

noCtrPts   = noPtsX * noPtsY * noPtsZ;
noDofs     = noCtrPts * 3;

% Boundary nodes
% 1. Displacement conditions:
%    - On bottom edge (y=0) with u_y = 0
%    - On right  edge (x=0) with u_x = 0

% 2. Traction conditions
%    - On left edge (x=-L) with t = (sigma_x,sigma_xy)
%    - On top edge  (y=L) with t = (sigma_xy,sigma_y)

% find boundary nodes for bounjdary conditions

leftNodes   = find(controlPts(:,1)==0)';
rightNodes  = find(controlPts(:,1)==a)';

% essential boundary conditions

uFixed     = zeros(size(leftNodes));
vFixed     = zeros(size(leftNodes));
wFixed     = zeros(size(leftNodes));

% build connectivity ...

generateIGA3DMesh

% initialization

K = sparse(noDofs,noDofs);  % global stiffness matrix
u = zeros(noDofs,1);        % displacement vector
f = zeros(noDofs,1);        % external force vector


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gauss quadrature rule

[W,Q] = gaussianQuadNURBS(p+1,q+1,r+1);

% Assembling system of equation
% Stiffness matrix and external force vector

disp([num2str(toc),'  ASSEMBLING THE SYSTEM'])

% Loop over elements (knot spans)

for e=1:noElems
    idu    = index(e,1);
    idv    = index(e,2);
    idw    = index(e,3);
    
    xiE    = elRangeU(idu,:); % [xi_i,xi_i+1]
    etaE   = elRangeV(idv,:); % [eta_j,eta_j+1]
    zetaE  = elRangeW(idw,:); % [zeta_k,zeta_k+1]
 
    sctr   = element(e,:);          %  element scatter vector
    nn     = length(sctr);
    
    sctrB(1:3:3*nn)    = 3*sctr-2;
    sctrB(2:3:3*nn)    = 3*sctr-1;
    sctrB(3:3:3*nn)    = 3*sctr-0;
    
    % loop over Gauss points
    
    for gp=1:size(W,1)
        pt      = Q(gp,:);
        wt      = W(gp);
        
        % compute coords in parameter space
        Xi      = parent2ParametricSpace(xiE,pt(1)  );
        Eta     = parent2ParametricSpace(etaE,pt(2) );
        Zeta    = parent2ParametricSpace(zetaE,pt(3));
        J2      = jacobianPaPaMapping3d(xiE,etaE,zetaE);
        
        % compute derivative of basis functions w.r.t parameter coord
        
        [N dRdxi dRdeta dRdzeta] = NURBS3DBasisDers([Xi;Eta;Zeta],...
                                   p,q,r,uKnot,vKnot,wKnot,weights');
        
        % compute the jacobian of physical and parameter domain mapping
        % then the derivative w.r.t spatial physical coordinates
        
        pts = controlPts(sctr,:);
        
        % Jacobian matrix
             
        jacob = pts'*[dRdxi' dRdeta' dRdzeta'];
        J1    = det(jacob);
        
        % Jacobian inverse and spatial derivatives
        
        invJacob   = inv(jacob);
        dRdx       = [dRdxi' dRdeta' dRdzeta'] * invJacob;
        
        % B matrix
        
        B = strainDispMatrix3d(nn,dRdx);
        
        % compute elementary stiffness matrix and
        % assemble it to the global matrix
        
        K(sctrB,sctrB) = K(sctrB,sctrB) + B' * D * B * J1 * J2 * wt;
    end
end

% Computing external force

f(3*21)=-10000;

disp([num2str(toc),'  APPLYING BOUNDARY  CONDITIONS'])

udofs = 3*leftNodes-2;    % global indecies  of the fixed x disps
vdofs = 3*leftNodes-1;    % global indecies  of the fixed y disps
wdofs = 3*leftNodes;      % global indecies  of the fixed z disps

[K,f]=applyDirichletBCs3D(K,f,udofs,vdofs,wdofs,uFixed,vFixed,wFixed);

% SOLVE SYSTEM

disp([num2str(toc),'  SOLVING THE SYSTEM'])
U=K\f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% POST-PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp([num2str(toc),'  POST-PROCESSING'])

Ux    = U(1:3:noDofs);
Uy    = U(2:3:noDofs);
Uz    = U(3:3:noDofs);

vtsFile = '../results/beam3D';
plotStress3d

figure
fac=5;
view(3)
hold on
plot_mesh(node+fac*[dispX dispY dispZ],elementV,'B8','r.-',1.2);


% figure 
% hold on
% plot_mesh(node,elementV,'B8','g.-');
% %plot_field(node,elementV,'B8',Uz); do not support B8 elements!!!
% view(3)
% plot3(controlPts(:,1),controlPts(:,2),controlPts(:,3),'r*');
% 
% buildVisualization3dMesh
% sigmaXX = zeros(size(node,1),1);
% sigmaYY = zeros(size(node,1),1);
% sigmaXY = zeros(size(node,1),1);
% 
% VTKPostProcess3d(node,elementV,'B8',...
%              [sigmaXX sigmaYY sigmaXY sigmaXX sigmaXX sigmaXX],[Ux Uy Uz]);



