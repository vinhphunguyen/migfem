%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Isogeometric analysis for two dimensional elasticity problems.
% 
% Curved beam in shear problem
% 
% Vinh Phu Nguyen, 
% Johns Hopkins University
% nvinhphu@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('~/code/xfem-efg-matlab/fem_util');
addpath C_files/

clc
clear all

global p q

refineCount   = 3; % 0: no refinement. Refine mesh with 1, 2 and so on 
noGPs         = 4; % # of Gauss points along one direction


E1  = 1000;  % Young modulus
E2  = 500;   % Young modulus
nu1 = 0.3;  % Poisson ratio
nu2 = 0.3;  % Poisson ratio
stressState ='PLANE_STRESS'; % either 'PLANE_STRAIN' or "PLANE_STRESS
F = 10;

% COMPUTE ELASTICITY MATRIX
C1 = elasticityMatrix(E1,nu1,stressState);
C2 = elasticityMatrix(E2,nu2,stressState);

tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%  CAD input: control points and knot vectors
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

annularBiMaterialData

% h-refinement here

if (refineCount) 
    hRefinement2d 
end

noCtrPts       = noPtsX * noPtsY;
noDofs         = noCtrPts * 2;

% find boundary nodes for bounjdary conditions

fixedXNodes  = [find(controlPts(:,2)==0);...
                find(controlPts(:,1)==0)];
fixedYNodes  =  find(controlPts(:,2)==0);
rightNodes   =  find(controlPts(:,1)==0)';

% plot the mesh

plotMesh (controlPts,weights,uKnot,vKnot,p,q,50,'r--','try.eps');

% build connectivity ...

generateIGA2DMesh

rightPoints    = controlPts(rightNodes,:);
rightEdgeMesh  = zeros(noElemsV,q+1);

for i=1:noElemsV
    if (q==1)    
      rightEdgeMesh(i,:) = [rightNodes(i) rightNodes(i+1)];
    else
      rightEdgeMesh(i,:) = [rightNodes(i) rightNodes(i+1) rightNodes(i+2)];
    end
end

% initialization

K = sparse(noDofs,noDofs);  % global stiffness matrix 
u = zeros(noDofs,1);        % displacement vector
f = zeros(noDofs,1);        % external force vector

jacob = zeros(2,2);

% essential boundary conditions

uFixed     = zeros(size(fixedXNodes))';
vFixed     = zeros(size(fixedYNodes))';

udofs      = fixedXNodes;          
vdofs      = fixedYNodes+noCtrPts;  

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
   connU  = elConnU(idu,:);
   connV  = elConnV(idv,:);
   
   noFnsU = length(connU);
   noFnsV = length(connV);
   
   sctr   = element(e,:);          %  element scatter vector
   sctrB  = [sctr sctr+noCtrPts]; %  vector that scatters a B matrix
   nn     = length(sctr); 
   
   B      = zeros(3,2*nn);
 
   if ( xiE(2) <= 0.5 ) 
       C = C1;
   else
       C = C2;
   end
   
   % loop over Gauss points 
   
    for gp=1:size(W,1)                        
      pt      = Q(gp,:);                          
      wt      = W(gp);                            
      
      % compute coords in parameter space
      Xi      = parent2ParametricSpace(xiE,pt(1)); 
      Eta     = parent2ParametricSpace(etaE,pt(2)); 
      J2      = jacobianPaPaMapping(xiE,etaE);
      
      % derivatives of shape functions
      
      [dRdxi dRdeta] = NURBS2Dders([Xi; Eta],p,q,uKnot,vKnot,weights');

      % compute the jacobian of physical and parameter domain mapping
      % then the derivative w.r.t spatial physical coordinates
            
      pts = controlPts(sctr,:);
      
      % Jacobian matrix
      
      jacob(1,1) = dRdxi  * pts(:,1);
      jacob(1,2) = dRdeta * pts(:,1);
      jacob(2,1) = dRdxi  * pts(:,2);
      jacob(2,2) = dRdeta * pts(:,2);
      
      J1         = det(jacob);
      
      % Jacobian inverse and spatial derivatives
      
      invJacob   = inv(jacob);
      dRdx       = [dRdxi' dRdeta'] * invJacob;
      
      % B matrix
      %        _                                      _
      %        |  N_1,x  N_2,x  ...      0      0  ... |
      %  B  =  |      0      0  ... N_1,y  N_2,y  ... |
      %        |  N_1,y  N_2,y  ... N_1,x  N_2,x  ... |
      %        -                                      -
      
      B(1,1:nn)       = dRdx(:,1)';
      B(2,nn+1:2*nn)  = dRdx(:,2)';
      B(3,1:nn)       = dRdx(:,2)';
      B(3,nn+1:2*nn)  = dRdx(:,1)';
      
      % compute elementary stiffness matrix and
      % assemble it to the global matrix
      
      K(sctrB,sctrB) = K(sctrB,sctrB) + B' * C * B * J1 * J2 * wt;
    end
end

% Computing external force

[W1,Q1] = quadrature(2, 'GAUSS', 1 ); 

% Loop over elements along left edge = noElemsV

for e=1:noElemsV
   xiE   = elRangeV(e,:); % [xi_i,xi_i+1]
   conn  = elConnV(e,:);
   noFns = length(conn);
   
   sctry = rightEdgeMesh(e,:)+noCtrPts;
   
   % loop over Gauss points 
    for gp=1:size(W1,1)                        
      xi      = Q1(gp,:);                          
      wt      = W1(gp);                            
      Xi      = 0.5 * ( ( xiE(2) - xiE(1) ) * xi + xiE(2) + xiE(1));
      J2      = 0.5 * ( xiE(2) - xiE(1) ); 
      
      N       = [];
      dNdxi   = [];
        
      % compute derivative of basis functions w.r.t parameter coord
      
      for in=1:noFns
       [Ni,dNi]  = NURBSbasis (conn(in),q,Xi,vKnot,weights);
       N         = [N Ni];
       dNdxi     = [dNdxi dNi];
      end
      
      % compute the jacobian of physical and parameter domain mapping
      % then the derivative w.r.t spatial physical coordinates
      
      jacob1   = dNdxi * rightPoints(conn,:);
      J1       = norm (jacob1);
         
      f(sctry) = f(sctry) + N' * F * J1 * J2 * wt;  
    end
end

disp([num2str(toc),'  APPLYING BOUNDARY  CONDITIONS'])
bcwt=mean(diag(K)); % a measure of the average  size of an element in K
                    % used to keep the  conditioning of the K matrix

applyBC


% SOLVE SYSTEM
disp([num2str(toc),'  SOLVING THE SYSTEM'])
U=K\f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% POST-PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ux    = U(1:noCtrPts);
Uy    = U(noCtrPts+1:noDofs);

scale = 10;
deformedControlPts = controlPts + scale * [Ux Uy];
figure
plotMesh (deformedControlPts,weights,uKnot,vKnot,p,q,50,'r-','deformed.eps');


uXNode1 = U(rightNodes(end)+noCtrPts)



% Compute stresses at Gauss points
% For visualization: triangulation GPs in physical coordinates
% Only works for linear elements!!!

[W,Q]=quadrature(  1, 'GAUSS', 2 );
stress    = zeros(3,noElems*2*2);
stressPts = []; % GPs in physical coordinates
ind       = 0;

for e=1:noElems
   idu    = index(e,1);
   idv    = index(e,2);
   xiE    = elRangeU(idu,:); % [xi_i,xi_i+1]
   etaE   = elRangeV(idv,:); % [eta_j,eta_j+1]
   connU  = elConnU(idu,:);
   connV  = elConnV(idv,:);
   
   noFnsU = length(connU);
   noFnsV = length(connV);
   
   sctr   = element(e,:);          %  element scatter vector
   sctrB  = [sctr sctr+noCtrPts]; %  vector that scatters a B matrix
   nn     = length(sctr); 
   
   B      = zeros(3,2*nn);
 
   % loop over Gauss points 
   
    for gp=1:size(W,1)                        
      pt      = Q(gp,:);                          
      wt      = W(gp);                            
      
      % compute coords in parameter space
      Xi      = parent2ParametricSpace(xiE,pt(1)); 
      Eta     = parent2ParametricSpace(etaE,pt(2)); 
      J2      = jacobianPaPaMapping(xiE,etaE);
      
      
      Nxi     = [];
      Neta    = [];
      dNdxi   = [];
      dNdeta  = [];
      dRdxi   = [];
      dRdeta  = [];
      N       = [];
        
      % compute derivative of basis functions w.r.t parameter coord
      
      for in=1:noFnsU
       [Ni,dNi]  = NURBSbasis (connU(in),p,Xi,uKnot,weights);
       Nxi       = [Nxi Ni];
       dNdxi     = [dNdxi dNi];
      end
      
      for in=1:noFnsV
       [Ni,dNi]  = NURBSbasis (connV(in),q,Eta,vKnot,weights);
       Neta      = [Neta Ni];
       dNdeta    = [dNdeta dNi];
      end
      
      % derivate of R=Nxi*Neta w.r.t xi and eta
      % this is derivative of shape functions in FEM
      
      for j=1:noFnsV
          for i=1:noFnsU
              dRdxi  = [dRdxi  dNdxi(i) * Neta(j)];
              dRdeta = [dRdeta Nxi(i)   * dNdeta(j)];
              N      = [N      Nxi(i)   * Neta(j)];
          end
      end
      
      % compute the jacobian of physical and parameter domain mapping
      % then the derivative w.r.t spatial physical coordinates
            
      pts = controlPts(sctr,:);
      
      % Jacobian matrix
      
      jacob(1,1) = dRdxi  * pts(:,1);jacob(1,2) = dRdeta * pts(:,1);
      jacob(2,1) = dRdxi  * pts(:,2);jacob(2,2) = dRdeta * pts(:,2);
      
      strPoint   = N * pts;
      stressPts  = [stressPts;strPoint];
     
      J1         = det(jacob);
      
      % Jacobian inverse and spatial derivatives
      
      invJacob   = inv(jacob);
      dRdx       = [dRdxi' dRdeta'] * invJacob;
      
      % B matrix
      
      B(1,1:nn)       = dRdx(:,1)';
      B(2,nn+1:2*nn)  = dRdx(:,2)';
      B(3,1:nn)       = dRdx(:,2)';
      B(3,nn+1:2*nn)  = dRdx(:,1)';
      
      strain          = B*U(sctrB);      
      ind             = ind + 1;
      stress(1:3,ind) = C*strain;
    end
end
% 
% hold on
% plotMesh (controlPts,weights,uKnot,vKnot,p,q,50,'b-','try.eps');

%plot(stressPts(:,1),stressPts(:,2),'r.');
% 

figure
hold on
tri = delaunay(stressPts(:,1),stressPts(:,2));
plot_field(stressPts,tri,'T3',stress(1,:));
axis('equal');
xlabel('X');
ylabel('Y');
title('Sigma XX');
set(gcf,'color','white');
colorbar('vert');







