%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Finite Cell Method for two dimensional elasticity problems.
%
% Infinite plate with a circular hole
%
% Vinh Phu Nguyen, March 2012
% nvinhphu@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath ../fem_util/
addpath ../C_files/
addpath ../data/
addpath ../meshing/
addpath ../post-processing/
addpath ../fem-functions/
addpath ../analytical-solutions/

clc
clear all

global p q W Q


E           = 1e5;  % Young modulus
nu          = 0.3;  % Poisson ratio
stressState = 'PLANE_STRESS';
L           = 4; % length of plate

vtuFile     = '../results/infinitePlateHoleFCM';

% Elasticity matrix

if ( strcmp(stressState,'PLANE_STRESS') )
    C=E/(1-nu^2)*[ 1      nu          0;
        nu     1          0 ;
        0     0  0.5*(1-nu) ];
else
    C=E/(1+nu)/(1-2*nu)*[ 1-nu  nu     0;
        nu    1-nu   0;
        0     0  0.5-nu ];
end

tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%  CAD input: control points and knot vectors
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plateHoleData

noGPs          = p+1; % # of GPs, assume p>=q

% Boundary nodes
% 1. Displacement conditions (symmetry):
%    - On bottom edge (y=0) with u_y = 0
%    - On left  edge (x=0) with u_x = 0

% 2. Traction conditions
%    - On right edge (x=L) with t = (sigma_x,sigma_xy)
%    - On top edge  (y=L) with t = (sigma_xy,sigma_y)

% find boundary nodes for boundary conditions

bottomNodes = find(controlPts(:,2)==0)';
rightNodes  = find(controlPts(:,1)==L)';
leftNodes   = find(controlPts(:,1)==0)';
topNodes    = find(controlPts(:,2)==L)';

inactiveDofs = [];

for ie=1:length(inactiveElems)    
    e    = inactiveElems(ie);
    sctr = element(e,:);
    inactiveDofs = [inactiveDofs sctr];
end

inactiveDofs = unique(inactiveDofs);

% essential boundary conditions

uNodes = [leftNodes   inactiveDofs];
vNodes = [bottomNodes inactiveDofs];

uNodes = unique(uNodes);
vNodes = unique(vNodes);

uFixed     = zeros(size(uNodes));
vFixed     = zeros(size(vNodes));

udofs = uNodes;           % global indecies  of the fixed x disps
vdofs = vNodes+noCtrPts;  % global indecies  of the fixed y disps

% build a 1D meshes for the right and top edges over which
% a traction is applied

rightPoints   = controlPts(rightNodes,:);
topPoints     = controlPts(topNodes,:);
leftEdgeMesh  = zeros(noElemsV,q+1);
topEdgeMesh   = zeros(noElemsU,p+1);

for i=1:noElemsV
    rightEdgeMesh(i,:) = rightNodes(i:i+q);
end

for i=1:noElemsU
    topEdgeMesh(i,:) = topNodes(i:i+p);
end

% initialization

K = sparse(noDofs,noDofs);  % global stiffness matrix
U = zeros(noDofs,1);        % displacement vector
f = zeros(noDofs,1);        % external force vector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Assembling system of equation
% Stiffness matrix and external force vector

coord=[-1 -1;1 -1;1 1;-1 1];

disp([num2str(toc),'  ASSEMBLING THE SYSTEM'])
%%
% Loop over elements (knot spans)

gps = [];

for e=1:noElems    
    if (ismember(e,inactiveElems)) 
        continue 
    end
        
    idu    = index(e,1);
    idv    = index(e,2);
    xiE    = elRangeU(idu,:); % [xi_i,xi_i+1]
    etaE   = elRangeV(idv,:); % [eta_j,eta_j+1]
    
    sctr   = element(e,:);         %  element scatter vector
    sctrB  = [sctr sctr+noCtrPts]; %  vector that scatters a B matrix
    nn     = length(sctr);
    B      = zeros(3,2*nn);
    pts    = controlPts(sctr,:);
    
    % -----------------------------------------------
    % Choose Gauss quadrature rules for elements
    
    if (ismember(e,splitElems))     % split element
        W = [];
        Q = [];
        [aa] = hierarchicalGaussQuad(noGPs,ls(elementV(e,:)),coord,0);
    else
        [W,Q] = quadrature(noGPs,'GAUSS',2);
    end
    
    % loop over Gauss points
    
    for gp=1:size(W,1)
        pt      = Q(gp,:);
        wt      = W(gp);
        
        % compute coords in parameter space
        Xi      = parent2ParametricSpace(xiE, pt(1));
        Eta     = parent2ParametricSpace(etaE,pt(2));
        J2      = jacobianPaPaMapping(xiE,etaE);
        
        % compute derivative of basis functions w.r.t parameter coord
        
        [N dRdxi dRdeta] = NURBS2DBasisDers([Xi; Eta],p,q,uKnot,vKnot,weights');
        
        % compute the jacobian of physical and parameter domain mapping
        % then the derivative w.r.t spatial physical coordinates
        % Jacobian matrix
        
        jacob     = pts'*[dRdxi' dRdeta'];
        J1        = det(jacob);
        invJacob  = inv(jacob);
        dRdx      = [dRdxi' dRdeta'] * invJacob;
        x         = N * pts; % global coord of GP
        gps       = [gps;x];
        
        if (sqrt(x(1)^2+x(2)^2)-r > 0)
            alpha = 1;
        else
            alpha = 0;
        end
        
        % B matrix
        
        B(1,1:nn)       = dRdx(:,1)';
        B(2,nn+1:2*nn)  = dRdx(:,2)';
        B(3,1:nn)       = dRdx(:,2)';
        B(3,nn+1:2*nn)  = dRdx(:,1)';
        
        % compute elementary stiffness matrix and
        % assemble it to the global matrix
        
        K(sctrB,sctrB) = K(sctrB,sctrB) + B' * C * B * alpha * J1 * J2 * wt;
    end
end

% Computing external force
%%
[W1,Q1] = quadrature(8, 'GAUSS', 1 );

% Loop over elements along right edge = noElemsV

for e=1:noElemsV
    xiE   = elRangeV(e,:); % [xi_i,xi_i+1]
    conn  = elConnV(e,:);
    pts   = rightPoints(conn,:);
    sctrx = rightEdgeMesh(e,:);
    sctry = sctrx + noCtrPts;
    
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
        
        % compute the exact stresses
        
        str = exact_plate_hole(x,r);
        tx  = str(1);
        ty  = str(3);
        fac = J1 * J2 * wt;
        f(sctrx) = f(sctrx) + N' * tx * fac;
        f(sctry) = f(sctry) + N' * ty * fac;
    end
end

% loop over top edge

for e=1:noElemsU
    xiE   = elRangeU(e,:); % [xi_i,xi_i+1]
    conn  = elConnU(e,:);
    pts   = topPoints(conn,:);
    sctrx = topEdgeMesh(e,:);
    sctry = sctrx + noCtrPts;
    
    % loop over Gauss points
    for gp=1:size(W1,1)
        xi      = Q1(gp,:);
        wt      = W1(gp);
        Xi      = 0.5 * ( ( xiE(2) - xiE(1) ) * xi + xiE(2) + xiE(1));
        J2      = 0.5 * ( xiE(2) - xiE(1) );
        
        [N dNdxi] = NURBS1DBasisDers(Xi,q,vKnot,weights);
        
        x        = N     * pts; % global coord of GP
        jacob1   = dNdxi * pts;
        J11      = norm (jacob1);
        
        % compute the exact stresses
        
        str = exact_plate_hole(x,r);
        tx  = str(3) ;
        ty  = str(2);
        fac = J11 * J2 * wt;
        f(sctrx) = f(sctrx) + N' * tx * fac;
        f(sctry) = f(sctry) + N' * ty * fac;
    end
end

%%
disp([num2str(toc),'  APPLYING BOUNDARY  CONDITIONS'])

applyBC

% SOLVE SYSTEM

disp([num2str(toc),'  SOLVING THE SYSTEM'])

U = K\f;
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% POST-PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ux    = U(1:noCtrPts);
Uy    = U(noCtrPts+1:noDofs);

%%

% Plot mesh and enriched nodes to check
figure
hold on
plot_mesh(node,elementV,'Q4','b-');
% plot the circle
theta = 0:0.01:pi/2;
xo = xc + r*cos(theta) ;
yo = yc + r*sin(theta) ;
plot(xo,yo,'k-','Linewidth',1.9);
% plot elements cut by the circle
plot_mesh(node,elementV(splitElems,:),'Q4','r-');
plot_mesh(node,elementV(inactiveElems,:),'Q4','c*-');
plot(gps(:,1),gps(:,2),'+');


%%%%%%%%%%%%%%%%%%%%%
%plot stress field

plotStressFCM

% compare exact and numerical displacements

for i = 1:size(node,1)
    x           = node(i,:) ;
    [exact_stress exact_disp] = exact_solution_hole(x,r,E,nu);
    ux_exact(i) = exact_disp(1);
    uy_exact(i) = exact_disp(2);
end
% ----------------------------------

% --------------------------------------------
% Plot both exact and numerical deformed shape
fac=600;
figure
hold on
h = plot(node(:,1)+fac*dispX,node(:,2)+fac*dispY,'rs');
set(h,'MarkerSize',7);
h = plot(node(:,1)+fac*ux_exact',node(:,2)+fac*uy_exact','b*');
set(h,'MarkerSize',7);
title('Exact and numerical deformed shape')
legend('XIGA','Exact')
axis equal




