%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Isogeometric analysis for two dimensional elasticity problems.
%
% Timoshenko beam
%
% Vinh Phu Nguyen,
% Delft University of Technology/Johns Hopkins University
% nvinhphu@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath ../fem_util/
addpath ../C_files/
addpath ../data/
addpath ../meshing/
addpath ../post-processing/
addpath ../fem-functions/
addpath ../nurbs-util/
addpath ../analytical-solutions/

clc
clear all

global element controlPts p q uKnot vKnot weights index elRangeU elRangeV
global Q W Q1 W1 rho C noDofs noCtrPts elConnU elConnV Ke0 Me0 damping noPtsX noPtsY

refineCount   = 0; % 0: no refinement. Refine mesh with 1, 2 and so on

E0  = 3e7;  % Young modulus
nu0 = 0.3;  % Poisson ratio

computeStress = 0; % do not compute stress for visualization

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

%timoshenkoBeamC0Data   % C0 similar to standard FEM
timoshenkoBeamC1Data % C2 B-spline (cubic elems)

noGPs  = 3;
noGPs1 = noGPs;

P = 1000;
I = (1/12)*D^3;

% h-refinement here

if (refineCount)
    hRefinement2d
end

noCtrPts       = noPtsX * noPtsY;
noDofs         = noCtrPts * 2;

% find boundary nodes for bounjdary conditions

fixedNodes  = find(controlPts(:,1)==0)';
forcedNodes = find(controlPts(:,1)==48)';

% plot the mesh

plotMesh (controlPts,weights,uKnot,vKnot,p,q,50,'b-','try.eps');

% build connectivity ...

generateIGA2DMesh

% build a 1D mesh for the right edge over which
% a traction is applied

bndPoints   = controlPts(forcedNodes,:);

bndMesh = zeros(noElemsV,q+1);

for i=1:noElemsV     
    bndMesh(i,:) = forcedNodes(i:i+q);
end

% initialization

K = sparse(noDofs,noDofs);  % global stiffness matrix
u = zeros(noDofs,1);        % displacement vector
f = zeros(noDofs,1);        % external force vector

jacob = zeros(2,2);

% essential boundary conditions, exact displacement is used
% as the left edge

uFixed     = zeros(size(fixedNodes));
vFixed     = zeros(size(fixedNodes));

% for i=1:length(fixedNodes)
%    id        = fixedNodes(i);
%    x         = controlPts(id,1);
%    y         = controlPts(id,2);
%    uFixed(i) = P*y/(6*E0*I)*((6*L-3*x)*x+(2+nu0)*(y^2-D^2/4));
%    vFixed(i) = -P/(6*E0*I)*(3*nu0*y^2*(L-x)+(4+5*nu0)*D^2*x/4+...
%                (3*L-x)*x^2);
% end


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
    nn     = length(sctr);
    
    sctrB(1,1:2:2*nn) = 2*sctr-1;    
    sctrB(1,2:2:2*nn) = 2*sctr  ;
    
    pts    = controlPts(sctr,:);
    
    B      = zeros(3,2*nn);
    
    % loop over Gauss points
    
    for gp=1:size(W,1)
        pt      = Q(gp,:);
        wt      = W(gp);
        
        [R,dRdx,J]  = getShapeGrads2D(pt,xiE,etaE,pts);                       
        B           = getBmatrix2D(B,dRdx); 
              
        % compute elementary stiffness matrix and
        % assemble it to the global matrix
        
        K(sctrB,sctrB) = K(sctrB,sctrB) + B' * C * B * J * wt;
    end
end

% Computing external force

[W1,Q1] = quadrature(noGPs1, 'GAUSS', 1 );

% Loop over elements along right edge = noElemsV

xR=[];
for e=1:noElemsV
    xiE   = elRangeV(e,:); % [xi_i,xi_i+1]
    conn  = elConnV(e,:);    
    sctry = 2*bndMesh(e,:);
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
        
        xR=[xR;x];        
        f(sctry) = f(sctry) + N' * ty * J1 * J2 * wt;
    end
end

disp([num2str(toc),'  APPLYING BOUNDARY  CONDITIONS'])

udofs=2*fixedNodes-1;          % global indecies  of the fixed x displacements
vdofs=2*fixedNodes;  % global indecies  of the fixed y displacements

K1 = K;

[K,f]=applyDirichletBCs(K,f,udofs,vdofs,uFixed,vFixed);

% SOLVE SYSTEM

disp([num2str(toc),'  SOLVING THE SYSTEM'])
U=K\f;

eu = 0.5*U'*K1*U
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% POST-PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find solution at point (L,0) and compare with exact solution
% uy_exact = -0.0089

Ux    = U(1:2:noDofs);
Uy    = U(2:2:noDofs);

ptId  = intersect (find(controlPts(:,1)==L),find(controlPts(:,2)==-D/2));
%ptId  = intersect (find(controlPts(:,1)==48),find(controlPts(:,2)==D/2));

x = L;
y = -D/2;

exactUy = -P/(6*E0*I)*(3*nu0*y^2*(L-x)+(4+5*nu0)*D^2*x/4+...
    (3*L-x)*x^2)

exactUx = P*y/(6*E0*I)*((6*L-3*x)*x+(2+nu0)*(y^2-D^2/4))

if (~isempty(ptId)) % there is 1 control pt at (L,0)
    %format long
    numericalUx = Ux(ptId)
    numericalUy = Uy(ptId)
else
    % must compute the displacement at point (L,0)
    % (L,0) => (xi,eta): nonlinear system of equations
    % xi,eta: use NURBS interpolation
    
    xi  = 1;
    eta = 1;
    
    projcoord = nurb2proj(noPtsX*noPtsY, [Ux Uy], weights);
    dim       = size(projcoord,2);
    tem       = SurfacePoint(noPtsX-1,p,uKnot,noPtsY-1,q,vKnot,...
        projcoord,dim,xi,eta);
    %format long
    numUy = tem(2)/tem(3)
    numUx = tem(1)/tem(3)
end

% plot the vertical displacemnt along the midline (y=0)
% both exact and numerical solutions.
%% Assume there are control points along this line!!!

% pei=P/(6*E0*I);
% pts = find(controlPts(:,2)==-D/2);
% y=-D/2;
% for i=1:length(pts)
%     x=controlPts(pts(i),1);
%     coordx(i)=x;
%     v_num(i)=Uy(pts(i));
%     v(i) =-pei*(3*nu0*y*y*(L-x)+(4+5*nu0)*D*D*x/4+(3*L-x)*x*x);
% end;
%  
% figure
% plot(coordx,v,'r-','LineWidth',1.1);
% hold on;
% plot(coordx,v_num,'bo','LineWidth',1);
% 
% xlabel('x')
% ylabel('Vertical displacement v at midplane line');
% legend('Analytical solution','NURBS');
%%

vtuFile = '../results/timoshenkoBeam';

plotStress1

clear disp; clear stress;

% +++++++++++++++++++++++++++++++++++++
%    COMPUTE THE ENERGY NORM
% +++++++++++++++++++++++++++++++++++++

disp([' COMPUTING THE ENERGY NORM ... '])

computeNorms

noDofs

    
    
    
    
