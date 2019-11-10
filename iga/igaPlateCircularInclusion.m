%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Isogeometric analysis for two dimensional elasticity problems.
%
% Infinite plate with a circular inclusion
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

refineCount = 2; % 0: no refinement. Refine mesh with 1, 2 and so on
noGPs       = 4; % # of Gauss points along one direction

E1  = 1000;  % Young modulus of matrix
E2  = 1;     % Young modulus of inclusion
nu1 = 0.3;  % Poisson ratio
nu2 = 0.3;  % Poisson ratio
stressState ='PLANE_STRAIN'; % either 'PLANE_STRAIN' or "PLANE_STRESS
t0  = -10;
L   = 4;
force = 0; % applied displacement used

% COMPUTE ELASTICITY MATRIX
C1 = elasticityMatrix(E1,nu1,stressState);
C2 = elasticityMatrix(E2,nu2,stressState);

tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%  CAD input: control points and knot vectors
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plateInclusionC1Data

% h-refinement here

if (refineCount)
    hRefinement2d
end

noCtrPts       = noPtsX * noPtsY;
noDofs         = noCtrPts * 2;

% Boundary nodes
% 1. Displacement conditions:
%    - On bottom edge (y=0) with u_y = 0
%    - On right  edge (x=0) with u_x = 0

% 2. Traction conditions
%    - On left edge (x=-L) with t = (sigma_x,sigma_xy)
%    - On top edge  (y=L) with t = (sigma_xy,sigma_y)

% find boundary nodes for bounjdary conditions

bottomNodes = find(controlPts(:,2)==0)';
rightNodes  = find(controlPts(:,1)==0)';
leftNodes   = find(controlPts(:,1)==-L)';
topNodes    = find(controlPts(:,2)==L)';

% essential boundary conditions
% value ubar=0
uFixed     = zeros(size(rightNodes));
vFixed     = zeros(size(bottomNodes));

udofs=rightNodes;          % global indecies  of the fixed x displacements
vdofs=bottomNodes+noCtrPts;  % global indecies  of the fixed y displacements

if (force==0)
    ubar   = -0.04;
    uFixed = [uFixed ubar*ones(1,length(leftNodes))];
    udofs  = [udofs leftNodes];
end

% plot the mesh

plotMesh (controlPts,weights,uKnot,vKnot,p,q,90,'r-','try.eps');

% build connectivity ...

generateIGA2DMesh

% build a 1D meshes for the left edge over which
% a traction is applied

leftPoints    = controlPts(leftNodes,:);
leftEdgeMesh  = zeros(noElemsV,q+1);

for i=1:noElemsU/2
    leftEdgeMesh(i,:) = leftNodes(i:i+q);
end

% initialization

K = sparse(noDofs,noDofs);  % global stiffness matrix
u = zeros(noDofs,1);        % displacement vector
f = zeros(noDofs,1);        % external force vector

%jacob = zeros(2,2);

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
    sctrB  = [sctr sctr+noCtrPts]; %  vector that scatters a B matrix
    nn     = length(sctr);
    
    if ( etaE(2) <= 0.3 )
        C = C2;
    else
        C = C1;
    end
    
    B      = zeros(3,2*nn);
    
    % loop over Gauss points
    
    for gp=1:size(W,1)
        pt      = Q(gp,:);
        wt      = W(gp);
        
        % compute coords in parameter space
        Xi      = parent2ParametricSpace(xiE, pt(1));
        Eta     = parent2ParametricSpace(etaE,pt(2));
        J2      = jacobianPaPaMapping(xiE,etaE);
        
        % compute derivative of basis functions w.r.t parameter coord
        
        [dRdxi dRdeta] = NURBS2Dders([Xi; Eta],p,q,uKnot,vKnot,weights');
        
        % compute the jacobian of physical and parameter domain mapping
        % then the derivative w.r.t spatial physical coordinates
        
        pts = controlPts(sctr,:);
        
        % Jacobian matrix
        
        %       jacob(1,1) = dRdxi * pts(:,1);jacob(1,2) = dRdeta * pts(:,1);
        %       jacob(2,1) = dRdxi * pts(:,2);jacob(2,2) = dRdeta * pts(:,2);
        
        jacob = pts'*[dRdxi' dRdeta'];
        J1    = det(jacob);
        
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

[W1,Q1] = quadrature(3, 'GAUSS', 1 );

% Loop over elements along left edge = noElemsV
if (force==1)
    for e=1:noElemsU/2
        xiE   = elRangeU(e,:); % [xi_i,xi_i+1]
        conn  = elConnU(e,:);
        noFns = length(conn);
        
        sctrx = leftEdgeMesh(e,:);
         
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
            
            jacob1   = dNdxi * leftPoints(conn,:);
            J1       = norm (jacob1);
            f(sctrx) = f(sctrx) + N' * t0 * J1 * J2 * wt;
        end
    end
end

disp([num2str(toc),'  APPLYING BOUNDARY  CONDITIONS'])

% use penalty method to enforce contraints u_a = u_b
% due to duplicate nodes at the top left corner

w    = 5000;
sctr = [leftNodes(end-1) leftNodes(end)];
penaltyStiffness = w*[1 -1;-1 1];
K(sctr,sctr) = K(sctr,sctr) + penaltyStiffness;

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

% controlPts=controlPts(1:66,:);
% weights=weights(1:66);
% scale = 1;
% deformedControlPts = controlPts + scale * [Ux Uy];
% figure
% plotMesh (deformedControlPts,weights,uKnot,vKnot,p,q,70,'r-','deformed.eps');
%plotMesh (controlPts,weights,uKnot,vKnot,p,q,70,'b--','deformed.eps');


% figure
% hold on
% tri = delaunay(controlPts(:,1),controlPts(:,2));
% plot_field(controlPts,tri,'T3',Ux);
% axis('equal');
% title('Displacement in x direction');
% %set(gcf,'color','white');
% colorbar('vert');


% plot stress field

vtuFile     = '../results/plateCircleInclusion';
plotStress1

stressComp=1;
figure
clf
plot_field(node,elementV,'Q4',stress(:,:,stressComp));
hold on
colorbar
title('Stress in x direction')
axis off
%plot_mesh(node,elementV,'Q4','g.-');

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,fileName,opts)












