%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Finite Cell Method for two dimensional elasticity problems.
%
% Infinite plate with a circular inclusion
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

global p q W Q levelMax


E           = 1e5;  % Young modulus
nu          = 0.3;  % Poisson ratio
stressState = 'PLANE_STRESS';
L           = 4; % length of plate
levelMax    = 4;

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
%

% find boundary nodes for boundary conditions

bottomNodes = find(controlPts(:,2)==0)';
rightNodes  = find(controlPts(:,1)==L)';
leftNodes   = find(controlPts(:,1)==0)';
topNodes    = find(controlPts(:,2)==L)';

inactiveDofs = [];

% for ie=1:length(inactiveElems)    
%     e    = inactiveElems(ie);
%     sctr = element(e,:);
%     inactiveDofs = [inactiveDofs sctr];
% end
% 
% inactiveDofs = unique(inactiveDofs);

% essential boundary conditions

uNodes = [leftNodes   inactiveDofs];
vNodes = [bottomNodes inactiveDofs];

uNodes = unique(uNodes);
vNodes = unique(vNodes);
vNodes = [vNodes topNodes];

uFixed     = zeros(size(uNodes));
vFixed     = [zeros(length(bottomNodes),1); 0.4*ones(length(topNodes),1)]';

udofs = uNodes;           % global indecies  of the fixed x disps
vdofs = vNodes+noCtrPts;  % global indecies  of the fixed y disps

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

    for gp=1:length(W)
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
            alpha = 0.001;
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
disp([num2str(toc),'  POST-PROCESSING'])

scale = 10;
deformedControlPts = controlPts + scale * [Ux Uy];
figure
plotMesh (deformedControlPts,weights,uKnot,vKnot,p,q,70,'r-','deformed.eps');

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






