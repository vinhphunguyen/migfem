%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Isogeometric analysis for two dimensional elasticity problems.
% 
% Square plate in uniaxial tension. Served as a verification
% of IGA code.
%
% This file is old. See other iga*.m files for a better implementation.
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

global p q

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem dependent data
% 1. Material properties
% 2. Geometry
% 3. Boundaries

%% 1. Material properties
E0          = 250;  % Young modulus
nu0         = 0.3;  % Poisson’s ratio
stressState ='PLANE_STRAIN'; % either 'PLANE_STRAIN' or "PLANE_STRESS
force       = 0;   % disp control, force=1: force control
ubar        = 0.4;

%% 2. Geometry 

plateC2Data    % in ../data folder

refineCount   = 3; % 0: no refinement. Refine mesh with 1, 2 and so on 
noGPs         = 2; % # of Gauss points surface integral along one direction
noGPs1        = 3; % # of GPs for line integral



% COMPUTE ELASTICITY MATRIX

C = elasticityMatrix(E0,nu0,stressState);

tic;



F = 100;


noCtrPts       = noPtsX * noPtsY;
noDofs         = noCtrPts * 2;

% find boundary nodes for bounjdary conditions

fixedNodes  = find(controlPts(:,1)==0)';
bottomNodes = find(controlPts(:,2)==0)';
forcedNodes = find(controlPts(:,1)==1)';

% plot the mesh

plotMesh (controlPts,weights,uKnot,vKnot,p,q,50,'b-','try.eps');

% build connectivity ...

generateIGA2DMesh

% build a 1D meshes for the right and top edges over which
% a traction is applied

rightPoints    = controlPts(forcedNodes,:);
rightEdgeMesh  = zeros(noElemsV,q+1);


for i=1:noElemsV
    rightEdgeMesh(i,:) = forcedNodes(i:i+q);
end

% initialization

K = sparse(noDofs,noDofs);  % global stiffness matrix 
u = zeros(noDofs,1);        % displacement vector
f = zeros(noDofs,1);        % external force vector

jacob   = zeros(2,2);
Nxi     = zeros(1,p+1);
Neta    = zeros(1,q+1);
dNdxi   = zeros(1,p+1);
dNdeta  = zeros(1,q+1);

% essential boundary conditions

uFixed     = zeros(size(fixedNodes));
vFixed     = zeros(1);

udofs=fixedNodes;          % global indecies  of the fixed x displacements
vdofs=intersect(fixedNodes,bottomNodes)+noCtrPts;  

if (force==0)
    uFixed = [uFixed ubar*ones(1,length(forcedNodes))];
    udofs = [udofs forcedNodes];
end


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
   
   sctr   = element(e,:);         %  element scatter vector
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
      
      dRdxi   = [];
      dRdeta  = [];
        
      % compute derivative of basis functions w.r.t parameter coord
      
      for in=1:noFnsU
       [Ni,dNi]  = NURBSbasis (connU(in),p,Xi,uKnot,weights);
       Nxi(in)    = Ni;
       dNdxi(in)  = dNi;
      end
      
      for in=1:noFnsV
       [Ni,dNi]  = NURBSbasis (connV(in),q,Eta,vKnot,weights);
       Neta(in)   = Ni;
       dNdeta(in) = dNi;
      end
      
      % derivate of R=Nxi*Neta w.r.t xi and eta
      % this is derivative of shape functions in FEM
      % the following code is for B-spline only!!!
      
      for j=1:noFnsV
          for i=1:noFnsU
              dRdxi  = [dRdxi  dNdxi(i) * Neta(j)];
              dRdeta = [dRdeta Nxi(i)   * dNdeta(j)];
          end
      end
      
      % for NURBS use the following
      
      %dNdxi 
      %dNdeta
      
      %[dRdxi dRdeta] = NURBS2Dders([Xi; Eta],p,q,uKnot,vKnot,weights'); 
      
      % compute the jacobian of physical and parameter domain mapping
      % then the derivative w.r.t spatial physical coordinates
            
      pts = controlPts(sctr,:);
      
      % Jacobian matrix
      
      jacob(1,1) = dRdxi  * pts(:,1);
      jacob(1,2) = dRdeta * pts(:,1);
      jacob(2,1) = dRdxi  * pts(:,2);
      jacob(2,2) = dRdeta * pts(:,2);
      
      J1         = det(jacob);
      %if (J1<0) J1 = -J1;end
      
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

if (force==1)
    [W1,Q1] = quadrature(noGPs1, 'GAUSS', 1 );
    
    % Loop over elements along left edge = noElemsV
    
    for e=1:noElemsV
        xiE   = elRangeV(e,:); % [xi_i,xi_i+1]
        conn  = elConnV(e,:);
        noFns = length(conn);
        
        sctrx = rightEdgeMesh(e,:);
        
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
            
            x        = N     * rightPoints(conn,:); % global coord of GP
            jacob1   = dNdxi * rightPoints(conn,:);
            J1       = norm (jacob1);
            f(sctrx) = f(sctrx) + N' * F * J1 * J2 * wt;
        end
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

scale = 0.8;
deformedControlPts = controlPts + scale * [Ux Uy];
figure
plotMesh (deformedControlPts,weights,uKnot,vKnot,...
           p,q,50,'b-','deformed.eps');

figure
hold on
tri = delaunay(controlPts(:,1),controlPts(:,2));
plot_field(controlPts,tri,'T3',Ux);
axis('equal');
title('Displacement in x direction');
%set(gcf,'color','white');
colorbar('vert');

figure
hold on
plot_field(controlPts,tri,'T3',Uy);
axis('equal');
title('Displacement in y direction');
%set(gcf,'color','white');
colorbar('vert');

projcoord = nurb2proj(noPtsX*noPtsY, Ux, weights);


dim=size(projcoord,2);


tem = SurfacePoint(noPtsX-1,p,uKnot,noPtsY-1,q,vKnot,...
    projcoord,dim,1,0);

tem(1)/tem(2)

vtuFile ='f.vtu';
plotStress1 





