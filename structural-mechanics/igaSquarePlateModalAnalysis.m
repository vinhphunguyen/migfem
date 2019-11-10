%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Isogeometric analysis for Kirchoff plate problems.
%
% Rotation-free thin plates. Fully clamped or simply supported 
% rectangular plates. Modal analysis.
%
% Vinh Phu Nguyen,
% Cardiff University
% nvinhphu@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath ../fem_util/
addpath ../C_files/
addpath ../data/
addpath ../meshing/
addpath ../post-processing/
addpath ../fem-functions/
addpath ../meshing/
addpath ../nurbs-util/
addpath ../nurbs-geopdes/inst/

clc
clear all

global p q

tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%  CAD input: control points and knot vectors
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plateData

% constitutive matrix

D  = E*t^3/(12*(1-nu^2));
C  = D*[1 nu 0;nu 1 0; 0 0 0.5*(1-nu)];

rho = 8000;

I0  = rho*t;
I2  = rho*t^3/12;

numberOfModes     = 12;
numberOfModesPlot = 12;

% find boundary nodes for boundary conditions

EPS = 1e-8;
bottomNodes  =  find(abs(controlPts(:,2))  <EPS);
topNodes     =  find(abs(controlPts(:,2)-b)<EPS);
leftNodes    =  find(abs(controlPts(:,1))  <EPS);
rightNodes   =  find(abs(controlPts(:,1)-a)<EPS);

fixedNodes  =  unique([bottomNodes;rightNodes;topNodes;leftNodes]);

if clamped
nextToBotNodes = noPtsX+2:2*noPtsX-1;
nextToRgtNodes = 2*noPtsX-1:noPtsX:noPtsX*(noPtsY-1)-1;
nextToTopNodes = noPtsX*(noPtsY-2)+2:noPtsX*(noPtsY-1)-1;
nextToLefNodes = noPtsX+2:noPtsX:noPtsX*(noPtsY-2)+2;

nextNodes      = unique([nextToBotNodes';nextToRgtNodes';...
                         nextToTopNodes';nextToLefNodes']);

fixedNodes     = [fixedNodes; nextNodes(:)];
end

plot(controlPts(fixedNodes,1),controlPts(fixedNodes,2),...
    'bs','MarkerEdgeColor','r','MarkerSize',14);

% build connectivity ...

generateIGA2DMesh

noCtrPts       = noPtsX   * noPtsY;
noDofs         = noCtrPts * 1;

% initialization

K = sparse(noDofs,noDofs);  % global stiffness matrix
M = sparse(noDofs,noDofs);  % global stiffness matrix

% essential boundary conditions

uFixed     = zeros(size(fixedNodes))';
udofs      = fixedNodes;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gauss quadrature rule
noGPs = p+1;
noGpEle = noGPs^2;
[W,Q]=quadrature(  noGPs, 'GAUSS', 2 ); % noGPs x noGPs point quadrature

% Assembling system of equation
% Stiffness matrix and external force vector

disp([num2str(toc),'  ASSEMBLING THE SYSTEM']);

% Loop over elements (knot spans)

for e=1:noElems
    idu    = index(e,1);
    idv    = index(e,2);
    xiE    = elRangeU(idu,:); % [xi_i,xi_i+1]
    etaE   = elRangeV(idv,:); % [eta_j,eta_j+1]
    
    sctr   = element(e,:);          %  element scatter vector
    nn     = length(sctr);
    pts    = controlPts(sctr,:);
    
    % loop over Gauss points
    
    for gp=1:size(W,1)
        pt      = Q(gp,:);
        wt      = W(gp);
        
        % compute coords in parameter space
        Xi      = parent2ParametricSpace(xiE,pt(1));
        Eta     = parent2ParametricSpace(etaE,pt(2));
        J2      = jacobianPaPaMapping(xiE,etaE);
        
        % shape functions, first and second derivatives w.r.t natural coords
        
        [R dRdxi dRdeta dR2dxi dR2det dR2dxe] = ...
            NURBS2DBasis2ndDers([Xi; Eta],p,q,uKnot,vKnot,weights');
        
        % compute the jacobian of physical and parameter domain mapping
        % then the derivative w.r.t spatial physical coordinates
        
        jacob  = [dRdxi; dRdeta]          * pts; % 2x2 matrix
        jacob2 = [dR2dxi; dR2det; dR2dxe] * pts; % 3x2 matrix
        
        J1    = det(jacob);
        
        dxdxi = jacob(1,1); dydxi = jacob(1,2);
        dxdet = jacob(2,1); dydet = jacob(2,2);
        
        j33   = [dxdxi^2     dydxi^2     2*dxdxi*dydxi;
            dxdet^2     dydet^2     2*dxdet*dydet;
            dxdxi*dxdet dydxi*dydet dxdxi*dydet+dxdet*dydxi];
        
        % Jacobian inverse and spatial 1st and 2nd derivatives
        
        invJacob   = inv(jacob);
        dRdx       = invJacob*[dRdxi;dRdeta];
        dR2dx      = inv(j33)*([dR2dxi; dR2det; dR2dxe]-jacob2*dRdx);
        
        % B matrix
        
        B          = dR2dx;
        B(3,:)     = B(3,:)*2;
        
        % compute elementary stiffness matrix and
        % assemble it to the global matrix
        
        wip = J1 * J2 * wt;
        K(sctr,sctr) = K(sctr,sctr) + B' * C * B * wip;
        M(sctr,sctr) = M(sctr,sctr) + (I0* R' * R + I2 * dRdx' * dRdx )* wip;
    end
end

activeDof=setdiff([1:noCtrPts]',[fixedNodes]);

[modeShape,freq]=eigs(K(activeDof,activeDof),M(activeDof,activeDof),...
    numberOfModes,0);

freq =diag(freq);
nfreq=freq*rho*t*a^4/D;
nfreq=nfreq.^(0.25);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% POST-PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






