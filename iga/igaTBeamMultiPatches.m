%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Isogeometric analysis for two dimensional elasticity problems.
% Illustration of a multi-patch IGA code.
% Timoshenko beam
%
% Vinh Phu Nguyen, February 2012
% Delft University of Technology
% nvinhphu@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('~/code/xfem-efg-matlab/fem_util');
addpath ../C_files/
addpath ../data/
addpath ../meshing/
addpath ../post-processing/
addpath ../fem-functions/
addpath ../analytical-solutions/
addpath ../nurbs-util/
addpath ~/code/igafem-nguyen/

clc
clear all

global p q 

E0  = 3e7;  % Young modulus
nu0 = 0.3;  % Poisson ratio

% Elasticity matrix

C   = E0/(1-nu0^2)*[  1      nu0          0;
    nu0        1          0;
    0        0  (1-nu0)/2  ];
tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%  CAD input: control points and knot vectors
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tBeamTwoPatchesC1Data

noGPs  = 3;
noGPs1 = noGPs;

P = 1000;
I = (1/12)*D^3;

noCtrPts       = max(max(chan));
noDofs         = noCtrPts * 2;

% find boundary nodes for boundary conditions

fixedNodes  = find(patch1.controlPts(:,1)==0)';
forcedNodes = find(patch2.controlPts(:,1)==48)';
bndPoints   = patch2.controlPts(forcedNodes,:);
forcedNodes = chan(:,end);

% build a 1D mesh for the right edge over which
% a traction is applied

bndMesh = zeros(patch2.noElemsV,q+1);

for i=1:patch2.noElemsV
    bndMesh(i,:) = forcedNodes(i:i+q);
end

% initialization

K = sparse(noDofs,noDofs);  % global stiffness matrix
u = zeros(noDofs,1);        % displacement vector
f = zeros(noDofs,1);        % external force vector

% essential boundary conditions, exact displacement is used
% as the left edge

uFixed     = zeros(size(fixedNodes));
vFixed     = zeros(size(fixedNodes));

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

% Loop over patches

for ip=1:noPatches
    index      = patches(ip).index;
    elRangeU   = patches(ip).elRangeU;
    elRangeV   = patches(ip).elRangeV;
    element    = patches(ip).element;    
    elementL   = patches(ip).elementLocal;    
    uKnot      = patches(ip).uKnot;
    vKnot      = patches(ip).vKnot;
    controlPts = patches(ip).controlPts;
    weights    = patches(ip).weights;
    noElems    = patches(ip).noElemsU * patches(ip).noElemsV;
    
    % Loop over elements (knot spans)
    
    for e=1:noElems
        idu    = index(e,1);
        idv    = index(e,2);
        xiE    = elRangeU(idu,:); % [xi_i,xi_i+1]
        etaE   = elRangeV(idv,:); % [eta_j,eta_j+1]
        
        sctr   = element(e,:);         %  global element scatter vector
        sctrL  = elementL(e,:);        %  local to the patch element scatter vector
        sctrB  = [sctr sctr+noCtrPts]; %  vector that scatters a B matrix
        nn     = length(sctr);
        pts    = controlPts(sctrL,:);
        
        B      = zeros(3,2*nn);
        
        % loop over Gauss points
        
        for gp=1:size(W,1)
            pt      = Q(gp,:);
            wt      = W(gp);
            
            % compute coords in parameter space
            Xi      = parent2ParametricSpace(xiE,pt(1));
            Eta     = parent2ParametricSpace(etaE,pt(2));
            J2      = jacobianPaPaMapping(xiE,etaE);
            
            % compute derivatives of shape functions
            [R dRdxi dRdeta] = NURBS2DBasisDers([Xi; Eta],p,q,uKnot,vKnot,weights');
            
            % compute the jacobian of physical and parameter domain mapping
            % then the derivative w.r.t spatial physical coordinates
            
            jacob      = pts' * [dRdxi' dRdeta'];
            J1         = det(jacob);                                    
            invJacob   = inv(jacob);
            dRdx       = [dRdxi' dRdeta'] * invJacob;
            
            % B matrix

            B(1,1:nn)       = dRdx(:,1)';
            B(2,nn+1:2*nn)  = dRdx(:,2)';
            B(3,1:nn)       = dRdx(:,2)';
            B(3,nn+1:2*nn)  = dRdx(:,1)';
            
            % elementary stiffness matrix and assemble it to the global matrix
            
            K(sctrB,sctrB) = K(sctrB,sctrB) + B' * C * B * J1 * J2 * wt;
        end % end of loop on Gauss points
    end     % end of loop on elements  
end         % end of loop on patches

% Computing the external force vector
%%

[W1,Q1] = quadrature(noGPs1, 'GAUSS', 1 );

% Loop over elements along right edge = noElemsV

for e=1:patch2.noElemsV
    xiE   = patch2.elRangeV(e,:); % [xi_i,xi_i+1]
    conn  = patch2.elConnV(e,:);
    sctry = bndMesh(e,:) + noCtrPts;
    pts   = bndPoints(conn,:);
    
    % loop over Gauss points
    
    for gp=1:size(W1,1)
        xi      = Q1(gp,:);
        wt      = W1(gp);
        Xi      = 0.5 * ( ( xiE(2) - xiE(1) ) * xi + xiE(2) + xiE(1));
        J2      = 0.5 * ( xiE(2) - xiE(1) );
        
        [N dNdxi] = NURBS1DBasisDers(Xi,q,vKnot,weights);
        
        % compute the jacobian of physical and parameter domain mapping
        % then the derivative w.r.t spatial physical coordinates
        
        x        = N    *pts; % global coord of GP
        jacob1   = dNdxi*pts;
        J1       = norm (jacob1);
        ty       = -(P/(2*I))*((D*D)/4-x(1,2)^2);
        
        f(sctry) = f(sctry) + N' * ty * J1 * J2 * wt;
    end
end

disp([num2str(toc),'  APPLYING BOUNDARY  CONDITIONS'])

udofs=fixedNodes;          % global indecies  of the fixed x displacements
vdofs=fixedNodes+noCtrPts; % global indecies  of the fixed y displacements

[K,f]=applyDirichletBCs(K,f,udofs,vdofs,uFixed,vFixed);

% SOLVE SYSTEM

disp([num2str(toc),'  SOLVING THE SYSTEM'])
U=K\f;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% POST-PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find solution at point (L,0) and compare with exact solution
% uy_exact = -0.0089

Ux    = U(1:noCtrPts);
Uy    = U(noCtrPts+1:noDofs);

x     = L;
y     = D/2;

exactUy     = -P/(6*E0*I)*(3*nu0*y^2*(L-x)+(4+5*nu0)*D^2*x/4+(3*L-x)*x^2)
exactUx     = P*y/(6*E0*I)*((6*L-3*x)*x+(2+nu0)*(y^2-D^2/4))
numericalUx = Ux(noCtrPts)
numericalUy = Uy(noCtrPts)


%%
% Stress computation

vtuFile0 = '../results/tBeamMultiPatch';

% Loop over patches

for ip=1:noPatches
    index      = patches(ip).index;
    elRangeU   = patches(ip).elRangeU;
    elRangeV   = patches(ip).elRangeV;
    element    = patches(ip).element;    
    elementL   = patches(ip).elementLocal;    
    uKnot      = patches(ip).uKnot;
    vKnot      = patches(ip).vKnot;
    controlPts = patches(ip).controlPts;
    weights    = patches(ip).weights;
    noPtsX     = patches(ip).noPtsX;
    noPtsY     = patches(ip).noPtsY;
    noElems    = patches(ip).noElemsU * patches(ip).noElemsV;
    
    vtuFile    = strcat(vtuFile0,num2str(ip));
    
    plotStressMP
end








