%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Isogeometric analysis for two dimensional elasticity problems.
%
% Timoshenko beam
%
% Dirichlet BCs: least-square method
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
timoshenkoBeamC1Data % C1 B-spline (quadratic elems)
%timoshenkoBeamC2Data % C2 B-spline (cubic elems)
%timoshenkoBeamC3Data % C3 B-spline (quartic elems)

refineCount = 2;

if (refineCount)
    %hRefinement2d
    hRefinement2dUniform
    noCtrPts       = noPtsX * noPtsY;
    noDofs         = noCtrPts * 2;
end

noGPs  = 3;
noGPs1 = noGPs;

P = 1000;
I = (1/12)*D^3;

%%

% find boundary nodes for bounjdary conditions

eps = 1e-12;
xx  = controlPts(:,1);
yy  = controlPts(:,2);

leftNodes   = find(xx==0)';
bottomNodes = find(abs(yy+6)<=eps)';
rightNodes  = find(abs(xx-48)<=eps)';
topNodes    = find(abs(yy-6)<=eps)';

topNodes    = sort(topNodes, 'descend');
leftNodes   = sort(leftNodes,'descend');

dispNodes   = [bottomNodes rightNodes(2:end) ...
               topNodes(2:end) leftNodes(2:end-1)];
noDispNodes = length(dispNodes);


% plot the mesh

plotMesh (controlPts,weights,uKnot,vKnot,p,q,150,'b-','try.eps');

% build connectivity ...

generateIGA2DMesh

bndElement     = zeros((noElemsU+noElemsV)*2,q+1); % assume p=q!!!

bottomEdgeMesh  = zeros(noElemsU,p+1);
rightEdgeMesh   = zeros(noElemsV,q+1);
topEdgeMesh     = zeros(noElemsU,p+1);
leftEdgeMesh    = zeros(noElemsV,q+1);

for i=1:noElemsU         
    bottomEdgeMesh(i,:) = bottomNodes(i:i+p);
    topEdgeMesh(i,:)    = topNodes   (i:i+p);
end

for i=1:noElemsV         
    leftEdgeMesh(i,:)  = leftNodes (i:i+q);
    rightEdgeMesh(i,:) = rightNodes(i:i+q);
end

% build boundary mesh for assembly matrix A in Aq=b

for ie=1:noElemsU
    bndElement(ie,:) = [ie:ie+p];
end

start = bottomNodes(end);

for ie=1:noElemsV
    bndElement(ie+noElemsU,:) = [start:start+p];
    start = start + 1;
end

start = max(max(bndElement));

for ie=1:noElemsU
    bndElement(ie+noElemsU+noElemsV,:) = [start:start+p];
    start = start + 1;
end

start = max(max(bndElement));

for ie=1:noElemsV
    bndElement(ie+2*noElemsU+noElemsV,:) = [start:start+p];
    start = start + 1;
end

bndElement((noElemsU+noElemsV)*2,p+1)=1;

% essential boundary conditions, exact displacement is used
% on the boundary edges

disp([num2str(toc),'  LEAST SQUARE for DIRICHLET BCs'])

A  = zeros(noDispNodes,noDispNodes);
bx = zeros(noDispNodes,1);
by = zeros(noDispNodes,1);

noxC   = 8;

% Loop over boundary edges...
% bottom edge

for ie=1:noElemsU
    sctr   = bottomEdgeMesh(ie,:);
    pts    = controlPts(sctr,:);    
    sctrA  = bndElement(ie,:);
    xiE    = elRangeU(ie,:);
    xiArr  = linspace(xiE(1),xiE(2),noxC);
   
    for ic=1:noxC                
        xi = xiArr(ic);
        [N dNdxi] = NURBS1DBasisDers(xi,p,uKnot,weights);
        
        A(sctrA,sctrA) = A(sctrA,sctrA) + N'*N;
        
        pt        = N    *pts; 
        
        x = pt(1); y = pt(2);
        
        % exact displacements
        
        ux = P*y/(6*E0*I)*((6*L-3*x)*x+(2+nu0)*(y^2-D^2/4));
        uy = -P/(6*E0*I)*(3*nu0*y^2*(L-x)+(4+5*nu0)*D^2*x/4+...
            (3*L-x)*x^2);
        
        bx(sctrA) = bx(sctrA) + ux*N';
        by(sctrA) = by(sctrA) + uy*N';
    end
end

% right edge

for ie=1:noElemsV
    sctr   = rightEdgeMesh(ie,:);
    pts    = controlPts(sctr,:);    
    sctrA  = bndElement(ie+noElemsU,:);
    xiE    = elRangeV(ie,:);
    xiArr  = linspace(xiE(1),xiE(2),noxC);
   
    for ic=1:noxC                
        xi = xiArr(ic);
        [N dNdxi] = NURBS1DBasisDers(xi,q,vKnot,weights);
        
        A(sctrA,sctrA) = A(sctrA,sctrA) + N'*N;
        
        pt        = N    *pts; 
        
        x = pt(1); y = pt(2);
        
        % exact displacements
        
        ux = P*y/(6*E0*I)*((6*L-3*x)*x+(2+nu0)*(y^2-D^2/4));
        uy = -P/(6*E0*I)*(3*nu0*y^2*(L-x)+(4+5*nu0)*D^2*x/4+...
            (3*L-x)*x^2);
        
        bx(sctrA) = bx(sctrA) + ux*N';
        by(sctrA) = by(sctrA) + uy*N';
    end
end

% top edge

for ie=1:noElemsU
    sctr   = topEdgeMesh(ie,:);
    pts    = controlPts(sctr,:);    
    sctrA  = bndElement(ie+noElemsU+noElemsV,:);
    xiE    = elRangeU(ie,:);
    xiArr  = linspace(xiE(1),xiE(2),noxC);
   
    for ic=1:noxC                
        xi = xiArr(ic);
        [N dNdxi] = NURBS1DBasisDers(xi,p,uKnot,weights);
        
        A(sctrA,sctrA) = A(sctrA,sctrA) + N'*N;
        
        pt        = N    *pts; 
        
        x = pt(1); y = pt(2);
        
        % exact displacements
        
        ux = P*y/(6*E0*I)*((6*L-3*x)*x+(2+nu0)*(y^2-D^2/4));
        uy = -P/(6*E0*I)*(3*nu0*y^2*(L-x)+(4+5*nu0)*D^2*x/4+...
            (3*L-x)*x^2);
        
        bx(sctrA) = bx(sctrA) + ux*N';
        by(sctrA) = by(sctrA) + uy*N';
    end
end

% left edge

for ie=1:noElemsV
    sctr   = leftEdgeMesh(ie,:);
    pts    = controlPts(sctr,:);    
    sctrA  = bndElement(ie+noElemsU*2+noElemsV,:);
    xiE    = elRangeV(ie,:);
    xiArr  = linspace(xiE(1),xiE(2),noxC);
   
    for ic=1:noxC                
        xi = xiArr(ic);
        [N dNdxi] = NURBS1DBasisDers(xi,q,vKnot,weights);
        
        A(sctrA,sctrA) = A(sctrA,sctrA) + N'*N;
        
        pt        = N    *pts; 
        
        x = pt(1); y = pt(2);
        
        % exact displacements
        
        ux = P*y/(6*E0*I)*((6*L-3*x)*x+(2+nu0)*(y^2-D^2/4));
        uy = -P/(6*E0*I)*(3*nu0*y^2*(L-x)+(4+5*nu0)*D^2*x/4+...
            (3*L-x)*x^2);
        
        bx(sctrA) = bx(sctrA) + ux*N';
        by(sctrA) = by(sctrA) + uy*N';
    end
end

% solve the system Aq_x=bx and Aq_y=by

[LL UU] = lu(A);
qxTemp  = LL\bx;
qyTemp  = LL\by;
qx      = UU\qxTemp;
qy      = UU\qyTemp;

%%%%

uFixed     = qx';
vFixed     = qy';

% for i=1:length(leftNodes)
%    id        = leftNodes(i);
%    x         = controlPts(id,1);
%    y         = controlPts(id,2);
%    uFixed(i) = P*y/(6*E0*I)*((6*L-3*x)*x+(2+nu0)*(y^2-D^2/4));
%    vFixed(i) = -P/(6*E0*I)*(3*nu0*y^2*(L-x)+(4+5*nu0)*D^2*x/4+...
%                (3*L-x)*x^2);
% end

udofs=2*dispNodes-1;          % global indecies  of the fixed x displacements
vdofs=2*dispNodes;  % global indecies  of the fixed y displacements

%%

% initialization

K = sparse(noDofs,noDofs);  % global stiffness matrix
u = zeros(noDofs,1);        % displacement vector
f = zeros(noDofs,1);        % external force vector

%%

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
        B           = getBmatrix2D(dRdx); 
        
        % compute elementary stiffness matrix and
        % assemble it to the global matrix
        
        K(sctrB,sctrB) = K(sctrB,sctrB) + B' * C * B * J * wt;
    end
end

disp([num2str(toc),'  APPLYING BOUNDARY  CONDITIONS'])

K1 = K;

[K,f]=applyDirichletBCs(K,f,udofs,vdofs,uFixed,vFixed);

% SOLVE SYSTEM

disp([num2str(toc),'  SOLVING THE SYSTEM'])
U=K\f;

eu = 0.5*U'*K1*U
%%

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

clear disp; 

% +++++++++++++++++++++++++++++++++++++
%    COMPUTE THE ENERGY NORM
% +++++++++++++++++++++++++++++++++++++

disp([' COMPUTING THE ENERGY NORM ... '])

% you have to look at file computeNorms.m, line 70
% because that file works for both problems:
% Timoshenko beam and the plate with a circular hole problem.
computeNorms

noDofs

    
    
    
    
