%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Isogeometric analysis for two dimensional elasticity problems.
% 
% L-shaped specimen problem
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
addpath ../analytical-solutions/
addpath ../nurbs-util/

clc
clear all

global element controlPts p q uKnot vKnot weights index elRangeU elRangeV
global Q W Q1 W1 rho C noDofs noCtrPts elConnU elConnV Ke0 Me0 damping noPtsX noPtsY

refineCount   = 6; % 0: no refinement. Refine mesh with 1, 2 and so on 
noGPs         = 3; % # of Gauss points along one direction


E0  = 1;  % Young modulus
nu0 = 0.3;  % Poisson ratio
stressState ='PLANE_STRESS'; % either 'PLANE_STRAIN' or "PLANE_STRESS
F = 1;

% COMPUTE ELASTICITY MATRIX

C = elasticityMatrix(E0,nu0,stressState);

tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%  CAD input: control points and knot vectors
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%LShapedData         %C0 elements
LShapedC1Data        %C1 elements 

% h-refinement here

if (refineCount) 
    hRefinement2d 
end

noCtrPts       = noPtsX * noPtsY;
noDofs         = noCtrPts * 2;

% find boundary nodes for bounjdary conditions

fixedXNodes  = find(controlPts(:,1)==2*a)';
fixedYNodes  = find(controlPts(:,2)==2*a)';
leftNodes    = find(controlPts(:,1)==0)';

% plot the mesh

plotMesh (controlPts,weights,uKnot,vKnot,p,q,50,'r--','try.eps');

% build connectivity ...

generateIGA2DMesh

leftPoints    = controlPts(leftNodes,:);
leftEdgeMesh  = zeros(noElemsV,q+1);

for i=1:noElemsV
   leftEdgeMesh(i,:) = leftNodes(i:i+q);
end

% initialization

K = sparse(noDofs,noDofs);  % global stiffness matrix 
u = zeros(noDofs,1);        % displacement vector
f = zeros(noDofs,1);        % external force vector

% essential boundary conditions

uFixed     = zeros(size(fixedXNodes));
vFixed     = zeros(size(fixedYNodes));

udofs      = 2*fixedXNodes-1;          
vdofs      = 2*fixedYNodes;  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gauss quadrature rule
[W,Q]=quadrature(  noGPs, 'GAUSS', 2 ); % noGPs x noGPs point quadrature

% Assembling system of equation
% Stiffness matrix and external force vector

disp([num2str(toc),'  ASSEMBLING THE SYSTEM'])

% Loop over elements (knot spans)

for e=1:noElems
   idu    = index(e,1);
   idv    = index(e,2);
   xiE    = elRangeU(idu,:); % [xi_i,xi_i+1]
   etaE   = elRangeV(idv,:); % [eta_j,eta_j+1]
   
   sctr   = element(e,:);          %  element scatter vector   
   nn     = length(sctr); 
   
   sctrB(1,1:2:2*nn) = 2*sctr-1;    
   sctrB(1,2:2:2*nn) = 2*sctr  ;
   
   pts    = controlPts(sctr,:);
   
   B      = zeros(3,2*nn);
 
   % loop over Gauss points 
   
    for gp=1:size(W,1)                        
      pt      = Q(gp,:);                          
      wt      = W(gp);                            
      
      [R,dRdx,J]  = getShapeGrads2D(pt,xiE,etaE,pts);                       
      B           = getBmatrix2D(B,dRdx); 
      
      % compute elementary stiffness matrix and
      % assemble it to the global matrix
      
      K(sctrB,sctrB) = K(sctrB,sctrB) + B' * C * B * J * wt;
    end
end

% Computing external force

[W1,Q1] = quadrature(2, 'GAUSS', 1 ); 

% Loop over elements along left edge = noElemsV

for e=1:noElemsV
   xiE   = elRangeV(e,:); % [xi_i,xi_i+1]
   conn  = elConnV(e,:);   
   pts   = leftPoints(conn,:);
   sctrx = 2*leftEdgeMesh(e,:)-1;
   
   % loop over Gauss points 
    for gp=1:size(W1,1)                        
      xi      = Q1(gp,:);                          
      wt      = W1(gp);                            
      Xi      = 0.5 * ( ( xiE(2) - xiE(1) ) * xi + xiE(2) + xiE(1));
      J2      = 0.5 * ( xiE(2) - xiE(1) ); 
      
      [N dNdxi] = NURBS1DBasisDers(Xi,q,vKnot,weights);
      
      % compute the jacobian of physical and parameter domain mapping
      % then the derivative w.r.t spatial physical coordinates
      
      x        = N     * pts; % global coord of GP
      jacob1   = dNdxi * pts;
      J1       = norm (jacob1);
         
      f(sctrx) = f(sctrx) + N' * F * J1 * J2 * wt;  
    end
end

disp([num2str(toc),'  APPLYING BOUNDARY  CONDITIONS'])
bcwt=mean(diag(K)); % a measure of the average  size of an element in K
                    % used to keep the  conditioning of the K matrix

[K,f]=applyDirichletBCs(K,f,udofs,vdofs,uFixed,vFixed);


% SOLVE SYSTEM
disp([num2str(toc),'  SOLVING THE SYSTEM'])
U=K\f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% POST-PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ux    = U(1:2:noDofs);
Uy    = U(2:2:noDofs);

scale = 0.8;
deformedControlPts = controlPts + scale * [Ux Uy];
%figure
%plotMesh (deformedControlPts,weights,uKnot,vKnot,p,q,50,'r-','deformed.eps');


uXNode1 = U(1)

% plot stress

vtuFile = '../results/lShaped';

plotStress1









