%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Isogeometric analysis for three dimensional elasticity problems.
%
% Thick cylinder subjects to line load
% 1/2 model is analysed.
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

clc
clear all

global p q r

E0           = 3e6;  % Young modulus
nu0          = 0.3;  % Poisson’s ratio

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE COMPLIANCE MATRIX
D=zeros(6,6);
D(1:3,1:3)=E0/(1+nu0)/(1-2*nu0)*[ 1-nu0 nu0 nu0;
                                  nu0 1-nu0 nu0;
                                  nu0 nu0 1-nu0];
D(4:6,4:6)=E0/(1+nu0)*eye(3);

tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%  CAD input: control points and knot vectors
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
thickCylinderData
noGPs        = 4; % # of Gauss points along one direction

noCtrPts   = noPtsX * noPtsY * noPtsZ;
noDofs     = noCtrPts * 3;

% plot NURBS mesh

plotMesh3(controlPts,weights, uKnot,vKnot,wKnot,...
                   p,q,r,40,'r-','try.eps');

% Boundary nodes
% find boundary nodes for bounjdary conditions

xConsNodes  = find(controlPts(:,1)==0)';
yConsNodes2 = find(controlPts(:,2)==-R)';
yConsNodes1 = find(controlPts(:,2)==R)';
yConsNodes  = [yConsNodes2 yConsNodes1];

% essential boundary conditions

u0 = -1;

uFixed     = zeros(size(xConsNodes));
vFixed     = [zeros(size(yConsNodes2)) ...
              u0*ones(size(yConsNodes1))];


udofs = xConsNodes;             % global indecies  of the fixed x disps
vdofs = yConsNodes+noCtrPts;    % global indecies  of the fixed y disps


aa         = find (controlPts(:,3)==300);
bb         = find (controlPts(:,2)==R);
cc         = find (controlPts(:,1)==0);
forcedNode = intersect(intersect(aa,bb),cc);

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
[W,Q]=quadrature(  noGPs, 'GAUSS', 3 ); 

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
    sctrB  = [sctr sctr+noCtrPts sctr+2*noCtrPts]; % scatters a B matrix
    nn     = length(sctr);
    
    B      = zeros(6,3*nn);
    
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

%f(forcedNode+noCtrPts)=-1/2;

disp([num2str(toc),'  APPLYING BOUNDARY  CONDITIONS'])

bcwt=mean(diag(K)); % a measure of the average  size of an element in K
% used to keep the  conditioning of the K matrix

f=f-K(:,udofs)*uFixed';  % modify the  force vector
f=f-K(:,vdofs)*vFixed';

f(udofs) = bcwt*uFixed;
f(vdofs) = bcwt*vFixed;

K(udofs,:)=0;  % zero out the rows and  columns of the K matrix
K(vdofs,:)=0;

K(:,udofs)=0;
K(:,vdofs)=0;

K(udofs,udofs)=bcwt*speye(length(udofs));  % put ones*bcwt on the diagonal
K(vdofs,vdofs)=bcwt*speye(length(vdofs));


% SOLVE SYSTEM

disp([num2str(toc),'  SOLVING THE SYSTEM'])
U=K\f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% POST-PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp([num2str(toc),'  POST-PROCESSING'])

Ux    = U(1:noCtrPts);
Uy    = U(noCtrPts+1:2*noCtrPts);
Uz    = U(2*noCtrPts+1:noDofs);

vtsFile = '../results/thickCylinder';
plotStress3d





