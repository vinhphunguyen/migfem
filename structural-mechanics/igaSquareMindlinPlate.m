%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Isogeometric analysis for Mindlin plate problems.
% High order NURBS eliminates shear locking nicely.
%
% Fully clamped or simply supported rectangular plates.
%
% Each control point has three dofs(deflection,rota1,rota2).
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

k   = 5/6;

D   = E*t^3/(12*(1-nu^2));
Cb  = D*[1 nu 0;nu 1 0; 0 0 0.5*(1-nu)];
Cs  = E*t*k/2/(1+nu)*[1 0;0 1];

% find boundary nodes for boundary conditions

EPS = 1e-8;
bottomNodes  =  find(abs(controlPts(:,2))  <EPS);
topNodes     =  find(abs(controlPts(:,2)-b)<EPS);
leftNodes    =  find(abs(controlPts(:,1))  <EPS);
rightNodes   =  find(abs(controlPts(:,1)-a)<EPS);

fixedNodes  =  unique([bottomNodes;rightNodes;topNodes;leftNodes]);

plot(controlPts(fixedNodes,1),controlPts(fixedNodes,2),...
    'bs','MarkerEdgeColor','r','MarkerSize',14);

% build connectivity ...

generateIGA2DMesh

noCtrPts       = noPtsX   * noPtsY;
noDofs         = noCtrPts * 3;

% initialization

K = sparse(noDofs,noDofs);  % global stiffness matrix
f = zeros(noDofs,1);        % external force vector

% essential boundary conditions

uFixed     = zeros(size(fixedNodes))';
udofs      = 3*fixedNodes-2;

vFixed = [];
wFixed = [];
vdofs  = []; % rotation 1
wdofs  = []; % rotation 2 

if clamped
    vFixed  = zeros(size(fixedNodes))';
    wFixed  = zeros(size(fixedNodes))';
    vdofs   = 3*fixedNodes-1;
    wdofs   = 3*fixedNodes;
end

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
            
    nn3    = nn*3;
    sctrw  = 3*sctr-2;
    sctrB  = zeros(1,nn3);
    
    sctrB(1:3:nn3) = 3*sctr-2; % deflection
    sctrB(2:3:nn3) = 3*sctr-1; % rotation 1
    sctrB(3:3:nn3) = 3*sctr;   % rotation 2
    
    Bb     = zeros(3,3*nn);
    Bs     = zeros(2,3*nn);
    
    % loop over Gauss points
    
    for gp=1:size(W,1)
        pt      = Q(gp,:);
        wt      = W(gp);
        
        % compute coords in parameter space
        Xi      = parent2ParametricSpace(xiE,pt(1));
        Eta     = parent2ParametricSpace(etaE,pt(2));
        J2      = jacobianPaPaMapping(xiE,etaE);
        
        % shape functions, first derivatives w.r.t natural coords
        
        [R dRdxi dRdeta ] = ...
            NURBS2DBasisDers([Xi; Eta],p,q,uKnot,vKnot,weights');
        
        % compute the jacobian of physical and parameter domain mapping
        % then the derivative w.r.t spatial physical coordinates
        
        jacob  = [dRdxi; dRdeta]          * pts; % 2x2 matrix               
        J1     = det(jacob);
                   
        % Jacobian inverse and spatial 1st derivatives
        
        invJacob   = inv(jacob);
        dRdx       = invJacob*[dRdxi;dRdeta];        
        
        % bending and shear B matrices
        
        Bb(1,2:3:3*nn) = dRdx(1,:);
        Bb(2,3:3:3*nn) = dRdx(2,:);
        Bb(3,2:3:3*nn) = dRdx(2,:);
        Bb(3,3:3:3*nn) = dRdx(1,:);
        
        Bs(1,2:3:3*nn) = -R;
        Bs(1,1:3:3*nn) = dRdx(1,:);
        Bs(2,1:3:3*nn) = dRdx(2,:);
        Bs(2,3:3:3*nn) = -R;
        
        % compute elementary stiffness matrix and
        % assemble it to the global matrix
        
        K(sctrB,sctrB) = K(sctrB,sctrB) + ...
                         (Bb' * Cb * Bb + Bs' * Cs * Bs) * J1 * J2 * wt;
                     (Bb' * Cb * Bb + Bs' * Cs * Bs) * J1 * J2 * wt
        f(sctrw)       = f(sctrw) + q0 * R' * J1 * J2 * wt;
    end
end

%%
disp([num2str(toc),'  APPLYING BOUNDARY  CONDITIONS']);

bcwt=mean(diag(K)); % a measure of the average  size of an element in K
% used to keep the  conditioning of the K matrix

% modify the force vector
f        = f-K(:,udofs)*uFixed';  
f        = f-K(:,vdofs)*vFixed';  
f        = f-K(:,wdofs)*wFixed';  
f(udofs) = bcwt*uFixed;
f(vdofs) = bcwt*vFixed;
f(wdofs) = bcwt*wFixed;

% modify the stiffness matrix
K(udofs,:)    = 0;  K(:,udofs)    = 0;
K(vdofs,:)    = 0;  K(:,vdofs)    = 0;
K(wdofs,:)    = 0;  K(:,wdofs)    = 0;

K(udofs,udofs)=bcwt*speye(length(udofs));  
K(vdofs,vdofs)=bcwt*speye(length(vdofs)); 
K(wdofs,wdofs)=bcwt*speye(length(wdofs)); 

%% SOLVE SYSTEM
disp([num2str(toc),'  SOLVING THE SYSTEM']);
U=K\f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% POST-PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Visualization using a Q4 visualization mesh

buildVisualizationMesh;
%plot_mesh(node,elementV,'Q4','g.-',1);

disp   = zeros(noElems,size(elementV,2));

for e=1:noElems
    idu    = index(e,1);
    idv    = index(e,2);
    xiE    = elRangeU(idu,:); % [xi_i,xi_i+1]
    etaE   = elRangeV(idv,:); % [eta_j,eta_j+1]
    
    sctr   = element(e,:);         %  element scatter vector
    sctrw  = 3*sctr-2;
    pts    = controlPts(sctr,:);
    
    uspan = FindSpan(noPtsX-1,p,xiE(1), uKnot);
    vspan = FindSpan(noPtsY-1,q,etaE(1),vKnot);
    
    % loop over Gauss points
    
    gp = 1;
    for iv=1:2
        if (iv==2)
            xiE = sort(xiE,'descend');
        end
        for iu=1:2
            Xi  = xiE(iu);
            Eta = etaE(iv);
            [N dRdxi dRdeta] = NURBS2DBasisDersSpecial([Xi; Eta],...
                p,q,uKnot,vKnot,weights',[uspan;vspan]);
            
            disp(e,gp)    = N * U(sctrw);
            
            gp = gp +1;
        end
    end
end

X = zeros(4,noElemsV);
Y = zeros(4,noElemsV);
Z = disp';

for i = 1:size(elementV,1)
    sctr   = elementV(i,:);
    X(:,i) = node(sctr,1);
    Y(:,i) = node(sctr,2);
end

figure
fill3(X,Y,Z,Z);
colorbar
title('deflection')
axis on

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,fileName,opts)


% Exact solution check

wbar = min(Z(:))*D*1000/q0/a^4
wext = 1.26532  % CCCC
%wext = 4.06235; % SSSS







