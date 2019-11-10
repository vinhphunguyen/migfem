%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Isogeometric analysis for two dimensional elasticity problems.
% 
% Cylinder subject to inner pressure. Only a quarter is modeled.
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
addpath ../meshing/
addpath ../nurbs-util/
addpath ../nurbs-geopdes/inst/

clc
clear all

global element controlPts p q uKnot vKnot weights index elRangeU elRangeV
global Q W Q1 W1 rho C noDofs noCtrPts elConnU elConnV Ke0 Me0 damping noPtsX noPtsY


E0          = 3e7;  % Young modulus
nu0         = 0.25;  % Poisson ratio
stressState ='PLANE_STRESS'; % either 'PLANE_STRAIN' or "PLANE_STRESS
pressure    = 3e4; % inner pressure 

% COMPUTE ELASTICITY MATRIX

C = elasticityMatrix(E0,nu0,stressState);

tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%  CAD input: control points and knot vectors
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

annularDataGeopdes

noCtrPts       = noPtsX * noPtsY;
noDofs         = noCtrPts * 2;

% find boundary nodes for boundary conditions

fixedXNodes  =  find(controlPts(:,1)==0);
fixedYNodes  =  find(controlPts(:,2)==0);
forcedNodes  =  1:noPtsX:noCtrPts;

% build connectivity ...

generateIGA2DMesh

% build boundary mesh for force vector computation

bndPoints      = controlPts(forcedNodes,:);
rightEdgeMesh  = zeros(noElemsV,q+1);

for i=1:noElemsV
    rightEdgeMesh(i,:) = forcedNodes(i:i+q);
end

% initialization

K = sparse(noDofs,noDofs);  % global stiffness matrix 
u = zeros(noDofs,1);        % displacement vector
f = zeros(noDofs,1);        % external force vector

% essential boundary conditions

uFixed     = zeros(size(fixedXNodes))';
vFixed     = zeros(size(fixedYNodes))';

udofs      = 2*fixedXNodes-1;          
vdofs      = 2*fixedYNodes;  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gauss quadrature rule
noGPs         = q+1; % # of Gauss points along one direction

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
   
   B      = zeros(3,2*nn);
   pts    = controlPts(sctr,:);
   
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

[W1,Q1] = quadrature(q+1, 'GAUSS', 1 ); 

% Loop over elements along left edge = noElemsV

for e=1:noElemsV
   xiE   = elRangeV(e,:); % [xi_i,xi_i+1]
   conn  = elConnV(e,:);
   noFns = length(conn);
   
   sctr  = rightEdgeMesh(e,:);
   sctrx = 2*sctr-1;
   sctry = 2*sctr;
   
   % loop over Gauss points 
    for gp=1:size(W1,1)                        
      xi      = Q1(gp,:);                          
      wt      = W1(gp);                            
      Xi      = 0.5 * ( ( xiE(2) - xiE(1) ) * xi + xiE(2) + xiE(1));
      J2      = 0.5 * ( xiE(2) - xiE(1) ); 
      
     [N dNdxi] = NURBS1DBasisDers(Xi,q,vKnot,weights);
      
      % compute the jacobian of physical and parameter domain mapping
      % then the derivative w.r.t spatial physical coordinates
      
      jacob1   = dNdxi * bndPoints(conn,:);
      J1       = norm (jacob1);
      
      x        = N *bndPoints(conn,:); % global coord of GP
      r        = norm(x);
      Fx       = pressure * x(1,1)/r;
      Fy       = pressure * x(1,2)/r;
      
      f(sctrx) = f(sctrx) + N' * Fx * J1 * J2 * wt;     
      f(sctry) = f(sctry) + N' * Fy * J1 * J2 * wt;  
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

vtuFile     = '../results/cylinderPressure';
plotStress1

stressComp=3;
figure
clf
plot_field(node,elementV,'Q4',stress(:,:,stressComp));
hold on
colorbar
title('Stress in x direction')
axis off
%plot_mesh(node,elementV,'Q4','g.-');

figure
clf
plot_field(node,elementV,'Q4',disp(:,:,1));
hold on
colorbar
title('Displacement in x direction')
axis off

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,fileName,opts)





