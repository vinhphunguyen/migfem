%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Isogeometric analysis for two dimensional elasticity problems.
%
% Timoshenko beam with boundary conditions imposed by the Lagrange
% multiplier method.
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
addpath ../analytical-solutions/
addpath ../nurbs-util/

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
%timoshenkoBeamC1Data % C1 B-spline (quadratic elems)
timoshenkoBeamC2Data % C2 B-spline (cubic elems)

noGPs  = 5;
noGPs1 = noGPs;

P = 1000;
I = (1/12)*D^3;

% find boundary nodes for bounjdary conditions

leftNodes   = find(controlPts(:,1)==0);
forcedNodes = find(controlPts(:,1)==48)';

% plot the mesh

plotMesh (controlPts,weights,uKnot,vKnot,p,q,50,'b-','try.eps');

% build connectivity ...

generateIGA2DMesh

% build a 1D mesh for the right edge over which
% a traction is applied

bndPoints   = controlPts(forcedNodes,:);
leftPoints  = controlPts(leftNodes,:);

bndMesh      = zeros(noElemsV,q+1);
leftEdgeMesh = zeros(length(leftNodes)-1,2);

for i=1:noElemsV     
    bndMesh(i,:) = forcedNodes(i:i+q);
end

for i=1:length(leftNodes)-1    
    leftEdgeMesh(i,:)  = leftNodes(i:i+1);    
end

% data initialization

K = sparse(noDofs,noDofs);  % global stiffness matrix
u = zeros(noDofs,1);        % displacement vector
f = zeros(noDofs,1);        % external force vector

disp_nodes     = leftNodes;
num_disp_nodes = length(disp_nodes);

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

% Computing external force

[W1,Q1] = quadrature(noGPs1, 'GAUSS', 1 );

% Loop over elements along right edge = noElemsV

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
           
        f(sctry) = f(sctry) + N' * ty * J1 * J2 * wt;
    end
end

%% Lagrange multiplier method to impose non-homogeneous bounadary conditions

disp([num2str(toc),'  APPLYING BOUNDARY  CONDITIONS'])

qk = zeros(1,2*num_disp_nodes);
G  = zeros(noDofs,2*num_disp_nodes);
    
% Loop over L2 finite elements (not NURBS elements) along 
% the left edge

m1 = 0 ;

for i = 1 : length(leftEdgeMesh)
    sctr = leftEdgeMesh(i,:);
    m1   = m1 + 1 ;
    m2   = m1 + 1 ;
    pts  = controlPts(sctr,:);  
    
    for gp = 1:size(W1,1)
        pt        = Q1(gp,:);
        wt        = W1(gp);
        [N,dNdxi] = lagrange_basis('L2',pt);       
        J0        = dNdxi'*pts;
        detJ      = norm(J0) ;
        pt        = N' * pts; % global GP
        
        % compute exact displacement
        
        x  = pt(1);
        y  = pt(2);
        
        ux = P*y/(6*E0*I)*((6*L-3*x)*x+(2+nu0)*(y^2-D^2/4));
        uy = -P/(6*E0*I)*(3*nu0*y^2*(L-x)+(4+5*nu0)*D^2*x/4+...
                (3*L-x)*x^2);
        
        N1 = N(1) ; N2 = N(2) ;
        
        % qk vector
        fac1       = wt * detJ * N1;
        fac2       = wt * detJ * N2;

        qk(2*m1-1) = qk(2*m1-1) - fac1 * ux;
        qk(2*m1)   = qk(2*m1)   - fac1 * uy;
        qk(2*m2-1) = qk(2*m2-1) - fac2 * ux;
        qk(2*m2)   = qk(2*m2)   - fac2 * uy;
        
        % G matrix
        
        xi        = inverseMapping1DNURBS (y,vKnot,leftPoints,q,2,elConnV,1);
        uspan     = FindSpan(noPtsY-1,q,xi,vKnot);
        nurbsCon  = leftNodes([uspan-q+1:uspan+1]);
        [N dNdxi] = NURBS1DBasisDers(xi,q,vKnot,weights);
                
        for j = 1 : length(nurbsCon)
            row1 = 2*nurbsCon(j)-1 ;
            row2 = 2*nurbsCon(j)   ;
            G1   = - wt * detJ * N(j) * [N1 0 ; 0 N1];
            G2   = - wt * detJ * N(j) * [N2 0 ; 0 N2];
            G(row1:row2,2*m1-1:2*m1) = G(row1:row2,2*m1-1:2*m1) + G1;
            G(row1:row2,2*m2-1:2*m2) = G(row1:row2,2*m2-1:2*m2) + G2;
        end
    end
end

% SOLVE SYSTEM
disp([num2str(toc),'  SOLVING THE SYSTEM'])

f = [f;qk'];                                % f = {f;qk}
m = ([K G; G' zeros(num_disp_nodes*2)]);    % m = [K GG;GG' 0]
d = m\f;                                    % d = {u;lamda}

% just get nodal parameters u_i, not need Lagrange multipliers
U = d(1:noDofs);
clear d;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% POST-PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find solution at point (L,0) and compare with exact solution
% uy_exact = -0.0089

Ux    = U(1:noCtrPts);
Uy    = U(noCtrPts+1:noDofs);

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

vtuFile = 'timoshenkoBeam';

plotStress1

clear disp; clear stress;

% +++++++++++++++++++++++++++++++++++++
%    COMPUTE THE ENERGY NORM
% +++++++++++++++++++++++++++++++++++++

disp([' COMPUTING THE ENERGY NORM ... '])

% you have to look at file computeNorms.m, line 70
% because that file works for both problems:
% Timoshenko beam and the plate with a circular hole problem.
computeNorms

noDofs

    
    
    
    
