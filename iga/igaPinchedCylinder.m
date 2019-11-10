%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Isogeometric analysis for three dimensional elasticity problems.
%
% Pinched cylinder subjects to point load
% 1/8 model is analysed.
%
% Vinh Phu Nguyen,
% Johns Hopkins University
% nvinhphu@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath ../fem_util/
addpath ../C_files/
addpath ../data/
addpath ../meshing/
addpath ../post-processing/
addpath ../fem-functions/
addpath ../nurbs-util/
addpath ../integration/

clc
clear all

global p q r uKnot vKnot wKnot

E0           = 3e6;  % Young modulus
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
    
%pinchedCylinderData % old data with simple refinement
pinchedCylinderCkData % new data using NURBS toolbox
                      % that allows k-refinement 

noGPs      = 4; % # of Gauss points along one direction

noCtrPts   = noPtsX * noPtsY * noPtsZ;
noDofs     = noCtrPts * 3;

%% Boundary nodes
% find boundary nodes for bounjdary conditions

rigidNodes  = find(controlPts(:,3)==0)'; 

% uncomment the following for free ends boundary conditions
%rigidNodes  = [];

zConsNodes  = find(abs(controlPts(:,3)-L) <= 1e-14)';
yConsNodes  = find(controlPts(:,2)==0)';
xConsNodes  = find(controlPts(:,1)==0)';

xConsNodes  = [xConsNodes rigidNodes];
yConsNodes  = [yConsNodes rigidNodes];
zConsNodes  = [zConsNodes];

xConsNodes  = unique(xConsNodes);
yConsNodes  = unique(yConsNodes);
zConsNodes  = unique(zConsNodes);

% essential boundary conditions

uFixed     = zeros(size(xConsNodes));
vFixed     = zeros(size(yConsNodes));
wFixed     = zeros(size(zConsNodes));

udofs      = 3*xConsNodes - 2;    % global indecies  of the fixed x disps
vdofs      = 3*yConsNodes - 1;    % global indecies  of the fixed y disps
wdofs      = 3*zConsNodes;        % global indecies  of the fixed z disps

% find loaded point index

aa         = find (controlPts(:,3)==L);
bb         = find (controlPts(:,2)==R);
cc         = find (controlPts(:,1)==0);
forcedNode = intersect(intersect(aa,bb),cc);
controlPts(forcedNode);
%% build connectivity ...

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
%[W,Q]=quadrature(  noGPs, 'GAUSS', 3 ); 
[W,Q]=gaussianQuadNURBS(p+1,q+1,r+1); 

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
    pts    = controlPts(sctr,:);
    B      = zeros(6,3*nn);
    
    % loop over Gauss points
    
    for gp=1:size(W,1)
        pt      = Q(gp,:);
        wt      = W(gp);
        
        % compute coords in parameter space
        Xi      = parent2ParametricSpace(xiE,  pt(1));
        Eta     = parent2ParametricSpace(etaE, pt(2));
        Zeta    = parent2ParametricSpace(zetaE,pt(3));
        J2      = jacobianPaPaMapping3d(xiE,etaE,zetaE);
        
        % compute derivative of basis functions w.r.t parameter coord
        
        [N dRdxi dRdeta dRdzeta] = NURBS3DBasisDers([Xi;Eta;Zeta],...
                                   p,q,r,uKnot,vKnot,wKnot,weights');
        
        % compute the jacobian of physical and parameter domain mapping
        % then the derivative w.r.t spatial physical coordinates                
        % Jacobian matrix
             
        jacob      = pts'*[dRdxi' dRdeta' dRdzeta'];
        J1         = det(jacob);                        
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

f(3*forcedNode-2)=-1/4;

disp([num2str(toc),'  APPLYING BOUNDARY  CONDITIONS'])

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

% vertical displacement at the loaded point
% The exact value is: 1.8248e-5 (rigid diaphrams)
% The exact value is: 4.52e-4   (free ends)

Uy(forcedNode)

vtsFile = '../results/pinchedCylinder';
plotStress3d





