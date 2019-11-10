%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Isogeometric analysis for two dimensional elasticity problems.
%
% Cylinder subject to inner pressure. Only a quarter is modeled.
%
% Fast assembly using the triplet sparse format.
%
% Vinh Phu Nguyen,
% Cardiff University, UK
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

tic;

%% Problem dependent data
%1. Material properties 
%2. Geometry 
%3. Force, output file name

%% 1. Material properties and force

E0          = 3e7;  % Young modulus
nu0         = 0.25;  % Poisson ratio
stressState ='PLANE_STRESS'; % either 'PLANE_STRAIN' or "PLANE_STRESS
pressure    = 3e4; % inner pressure
vtuFile     = '../results/cylinderPressure';

%% 2. Geometry

annularDataGeopdes   % annularDataGeopdes.m 


%% Processing ...

% initialization

u = zeros(noDofs,1);        % displacement vector
f = zeros(noDofs,1);        % external force vector

%% one the assembly
nElNod = size(element,2);
nElDof = nElNod*2;
nElmLK = nElDof^2;
nSprGK = nElmLK*noElems;

jSprLK = 1:nElmLK;
vIdGrd = ones(1,nElDof);

% values, row indices, columns indices of the global K matrix
vSprGK = zeros(nSprGK,1);
jSprRw = zeros(nSprGK,1);
jSprCl = zeros(nSprGK,1);

Ke0    = zeros(nElDof,nElDof); % element Ke

%% essential boundary conditions

uFixed     = zeros(size(fixedXNodes))';
vFixed     = zeros(size(fixedYNodes))';

udofs      = 2*fixedXNodes-1;
vdofs      = 2*fixedYNodes;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% COMPUTE ELASTICITY MATRIX

C = elasticityMatrix(E0,nu0,stressState);

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
    sctrB(1,2:2:2*nn) = 2*sctr;
    
    B      = zeros(3,2*nn);
    
    pts    = controlPts(sctr,:);
    
    Ke     = Ke0;
    
    % loop over Gauss points
    
    for gp=1:size(W,1)
        pt      = Q(gp,:);
        wt      = W(gp);
        
        % compute coords in parameter space
        Xi      = parent2ParametricSpace(xiE,pt(1));
        Eta     = parent2ParametricSpace(etaE,pt(2));
        J2      = jacobianPaPaMapping(xiE,etaE);
        
        [dRdxi dRdeta] = NURBS2Dders([Xi; Eta],p,q,uKnot,vKnot,weights');
        
        % compute the jacobian of physical and parameter domain mapping
        % then the derivative w.r.t spatial physical coordinates
        
        jacob = pts'*[dRdxi' dRdeta'];
        J1    = det(jacob);
        
        % Jacobian inverse and spatial derivatives
        
        invJacob   = inv(jacob);
        dRdx       = [dRdxi' dRdeta'] * invJacob;
        
        % B matrix
        B(1,1:2:2*nn)  = dRdx(:,1)';
        B(2,2:2:2*nn)  = dRdx(:,2)';
        B(3,1:2:2*nn)  = dRdx(:,2)';
        B(3,2:2:2*nn)  = dRdx(:,1)';
        
        % compute elementary stiffness matrix and
        % assemble it to the global matrix
        
        Ke = Ke + B' * C * B * J1 * J2 * wt;
    end
    
    %----------------------------------------------------------------------
    % Global Assembly (Sparse)
    %----------------------------------------------------------------------
    
    vElDof = sctrB';
    mRwGrd = vElDof(:,vIdGrd);
    mClGrd = mRwGrd';
    
    jSprRw(jSprLK) = mRwGrd(:); % put columns of mRwGrd one after each other
    jSprCl(jSprLK) = mClGrd(:); % in a column vector 
    vSprGK(jSprLK) = Ke(:);
    
    jSprLK         = jSprLK + nElmLK; % move to the next position
end

% Here comes the total stiffness matrix in one shot!!!

K = sparse(jSprRw,jSprCl,vSprGK,noDofs,noDofs);

%% Computing external force

[W1,Q1] = quadrature(q+1, 'GAUSS', 1 );

% Loop over elements along left edge = noElemsV

for e=1:size(rightEdgeMesh,1)
    xiE   = neumannMesh.range(e,:); % [xi_i,xi_i+1]            
    sctr  = neumannMesh.element(e,:);
    sctrx = 2*sctr-1;
    sctry = 2*sctr;
    pts   = controlPts(sctr,:);
    % loop over Gauss points
    for gp=1:size(W1,1)
        xi      = Q1(gp,:);
        wt      = W1(gp);
        Xi      = 0.5 * ( ( xiE(2) - xiE(1) ) * xi + xiE(2) + xiE(1));
        J2      = 0.5 * ( xiE(2) - xiE(1) );
        
        [N dNdxi] = NURBS1DBasisDers(Xi,neumannMesh.p,neumannMesh.knot,weights);
        
        % compute the jacobian of physical and parameter domain mapping
        % then the derivative w.r.t spatial physical coordinates
        
        jacob1   = dNdxi * pts;
        J1       = norm (jacob1);
        
        x        = N * pts; % global coord of GP
        r        = norm(x);
        Fx       = pressure * x(1,1)/r;
        Fy       = pressure * x(1,2)/r;
        
        f(sctrx) = f(sctrx) + N' * Fx * J1 * J2 * wt;
        f(sctry) = f(sctry) + N' * Fy * J1 * J2 * wt;
    end
end

disp([num2str(toc),'  APPLYING BOUNDARY  CONDITIONS'])

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





