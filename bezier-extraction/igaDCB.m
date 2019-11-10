%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Isogeometric analysis for two dimensional delamination problems.
% Using the so-called Bezier extraction operator.
%
% Double cantilever beam.
%
% Vinh Phu Nguyen,
% April, 2013
% Cardiff University, UK
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

global element C Cxi Cet weights controlPts W Q W1 Q1 shapes derivs noDofs
global noElems De noBasis Gc strength dummy iElements p q noElemsU noElemsV
global damage0 damage

E0          = 100;       % Young modulus
nu0         = 0.3;          % Poisson ratio
stressState ='PLANE_STRAIN'; % either 'PLANE_STRAIN' or "PLANE_STRESS
Gc          = 5.8;          % fracture energy
strength    = 28.5;            % tensile strength
dummy       = 1e6;



deltaf      = 0.5;
noSteps     = 20;

% COMPUTE ELASTICITY MATRIX

De = elasticityMatrix(E0,nu0,stressState);

tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%  CAD input: control points and knot vectors
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dcbCkData

% find boundary nodes for boundary conditions

eps=1e-10;

fixedNodes   =  find(abs(controlPts(:,2)) <eps);
forcedNodes  =  find(abs(controlPts(:,2)-w) <eps);

% build connectivity ...

generateIGA2DMesh

% initialization


u    = zeros(noDofs,1);        % displacement vector
fext = zeros(noDofs,1);        % external force vector

% essential boundary conditions

uFixed     = zeros(size(fixedNodes))';
vFixed     = zeros(size(fixedNodes))';

udofs      = 2*fixedNodes-1;
vdofs      = 2*fixedNodes;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gauss quadrature rule

noGPs         = p+1; % # of Gauss points along one direction

[W,Q]   = quadrature( noGPs, 'GAUSS', 2 ); % noGPs x noGPs point quadrature
[W1,Q1] = quadrature( noGPs, 'GAUSS', 1 ); % for interface elements

damage0  = zeros(noElemsU*(p+1),1);         % damage history variables
damage   = zeros(noElemsU*(p+1),1);         % damage history variables

%% Pre-compute Bernstein basis and derivatives for ONE Bezier element

noBasis = (p+1)*(q+1);
noGpEle = noGPs*noGPs;

shapes  = zeros(noGpEle,noBasis);
derivs  = zeros(noGpEle,noBasis,2);

for gp=1:size(W,1)
    [shapes(gp,:) derivs(gp,:,:)] = getShapeGradBernstein2D(p,q,Q(gp,1),Q(gp,2));
end


tol      = 1e-4;
iterMax  = 20;  % maximum number of Newton-Raphson iterations

for i = 1:noSteps
    
    disp ('=================================')
    disp (sprintf('   Load step %d',i));
    disp ('=================================')
    disp ('  NR iter : L2-norm residual')
    
    % compute the tangent and internal force
    
    [K,fint] = getStiffness2D (u);
    [K,fint] = getStiffnessForInterfaceElems (K,fint,u);
    
    % external force (the applied force is in vertical direction)
    
    %fext(2*forcedNodes1) = -deltaf;
    fext(2*forcedNodes) =  fext(2*forcedNodes) + deltaf;
    
    error    = 1;
    iiter    = 0;
    
    while error > tol
        
        iiter    = iiter + 1;
        
        % solve for the displacement increment
        
        res      = fext-fint;
        [K,res]  = applyDirichletBCs(K,res,udofs,vdofs,uFixed,vFixed);
        deltaU   = K\res;
        
        % update displacements
        
        u = u + deltaU;
        
        [K,fint] = getStiffness2D (u);
        [K,fint] = getStiffnessForInterfaceElems (K,fint,u);
        
        if iiter == 1
            rnmax = 0;
        end
        
        ndeltaU = norm(deltaU);
        rnmax = max (rnmax, ndeltaU);
        error =norm(deltaU);% ndeltaU/rnmax;
        
        disp (sprintf(' %s %i %s %5.4e ', 'Iter',iiter, ':', error) );
        
        if iiter == iterMax
            error('Newton-Raphson iterations did not converge!');
        end
    end
    
    %% the step has been converged...
    
    damage0 = damage;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% POST-PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ux    = u(1:2:noDofs);
Uy    = u(2:2:noDofs);

fac = 1e5;
plotMesh (controlPts+fac*[Ux Uy],weights,uKnot,vKnot,p,q,100,'r--','try.eps');


%% convert to NURBS object for visualization

defControlPts = controlPts+fac*[Ux Uy];

controlPts1          = zeros(4,noPtsX,noPtsY);
for i=1:noPtsY
    en  = noPtsX*(i-1)+1;
    st  = noPtsX*i;
    controlPts1(1:2,1:noPtsX,i) = defControlPts(en:st,:)';
end
controlPts1(4,:,:)   = 1;

solid = nrbmak(controlPts1,{uKnot vKnot});

nrbplot(solid,[80 80]);
view ([0 90])










