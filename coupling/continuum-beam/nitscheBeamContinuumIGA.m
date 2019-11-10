% This file implements the Nitsche method to join two mechanical models:
% a 2D continuum model and a Timoshenko beam model.
%
% Discretisation: B-spline elements.
%
% Problem: Timoshenko beam in bending.
%
% Vinh Phu Nguyen
% Cardiff University, UK
% 3 July 2013


addpath ../../gmshFiles/
addpath ../../post-processing/
addpath ../../fem-functions/
addpath ../../analytical-solutions/
addpath ../../fem_util/;
addpath ../../nurbs-geopdes/inst/
addpath ../../nurbs-util/
addpath ../../meshing/
addpath ../../fem-functions/
addpath ../../post-processing/
addpath ../../xiga/
addpath ../../nurbs-geopdes/inst/
addpath ../../delamination/
addpath ../../xml_toolbox/
addpath ../../C_files/


clear all
colordef white
state = 0;
tic;

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,'splinecurve.eps',opts)


% MATERIAL PROPERTIES
E0  = 3e7;  % Young?s modulus
nu0 = 0.3;  % Poisson?s ratio

% BEAM PROPERTIES
L  = 48;     % length of the beam
c  = 3;      % the distance of the outer fiber of the beam from the mid-line

t  = 2*c; % thickness, area = bxt
b  = 1; %
I0=2*c^3/3;  % the second polar moment of inertia of the beam cross-section.
P  = -1000*t^3/12/I0; % end force
G  = E0/2/(1+nu0);
k  = E0*I0;
shearFactor = 10*(1+nu0)/(12+11*nu0);

plotMesh  = 1;

stressState = 'PLANE_STRESS';
% COMPUTE ELASTICITY MATRIX
C = elasticityMatrix(E0,nu0,stressState);

% general averaged stress {sigma}=gamma1*sigma^s+gamma2*sigma^b
% with gamma1+gamma2=1

gamma2=0.5;
gamma1=1-gamma2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE FINITE ELEMENT MESH
%

disp([num2str(toc),'   GENERATING MESH'])

% MESH PROPERTIES

L1 = 24;
%% meshing for domain1-----------------------------------------------------

controlPts          = zeros(4,2,2);

controlPts(1:2,1,1) = [0;-c];
controlPts(1:2,2,1) = [L1;-c];
controlPts(1:2,1,2) = [0;c];
controlPts(1:2,2,2) = [L1;c];

controlPts(4,:,:)   = 1;

uKnot = [0 0 1 1];
vKnot = [0 0 1 1];

solid1 = nrbmak(controlPts,{uKnot vKnot});

% evaluate order
solid1    = nrbdegelev(solid1,[2 2]);

% knot insertion
uKnot     = cell2mat(solid1.knots(1));
vKnot     = cell2mat(solid1.knots(2));

refineCountX = 4;
for i=1:refineCountX
    uKnotVectorU = unique(uKnot);
    uKnotVectorV = unique(vKnot);    
    newKnotsX = uKnotVectorU(1:end-1) + 0.5*diff(uKnotVectorU);
    newKnotsY = uKnotVectorV(1:end-1) + 0.5*diff(uKnotVectorV);
    newKnots  = {[] newKnotsY};
    solid1     = nrbkntins(solid1,newKnots);
    uKnot      = cell2mat(solid1.knots(1));
    vKnot      = cell2mat(solid1.knots(2));
end

%solid1     = nrbkntins(solid1,{[] [0.2 0.4 0.6 0.8]});

refineCountX = 6;
for i=1:refineCountX
    uKnotVectorU = unique(uKnot);
    uKnotVectorV = unique(vKnot);
    newKnotsY = [];
    newKnotsX = uKnotVectorU(1:end-1) + 0.5*diff(uKnotVectorU);
    newKnotsY = uKnotVectorV(1:end-1) + 0.5*diff(uKnotVectorV);
    newKnots  = {newKnotsX []};
    solid1     = nrbkntins(solid1,newKnots);
    uKnot      = cell2mat(solid1.knots(1));
    vKnot      = cell2mat(solid1.knots(2));
end

mesh1 = buildIGA2DMesh (solid1);

%% domain 2----------------------------------------------------------------

controlPts        = zeros(4,2);
controlPts(1:2,1) = [L1;0];
controlPts(1:2,2) = [L;0];
controlPts(4,:)   = 1;

uKnot = [0 0 1 1];

solid2 = nrbmak(controlPts,uKnot);

% evaluate order
solid2 = nrbdegelev(solid2,2);

% knot insertion
uKnot     = solid2.knots;

refineCountX = 2;
for i=1:refineCountX
    uKnotVectorU = unique(uKnot);        
    newKnotsX  = uKnotVectorU(1:end-1) + 0.5*diff(uKnotVectorU);    
    newKnots   = newKnotsX;
    solid2     = nrbkntins(solid2,newKnots);
    uKnot      = solid2.knots;    
end

mesh2           = buildGA1DMesh (solid2.knots,solid2.order-1);
mesh2.elemCount = size(mesh2.element,1);
mesh2.controlPts= solid2.coefs(1:2,:)';

% PLOT MESH
if ( plotMesh )
    clf
    hold on
    nrbkntplot(solid1);
    plot(mesh1.controlPts(:,1),mesh1.controlPts(:,2),'o','MarkerEdgeColor','k',...
        'MarkerFaceColor','g','MarkerSize',6.5);
    nrbkntplot(solid2);
    plot(mesh2.controlPts(:,1),mesh2.controlPts(:,2),'o','MarkerEdgeColor','k',...
        'MarkerFaceColor','w','MarkerSize',6.5);
    hold on
    axis off
    axis([0 L -c c])
end

%% boundary mesh for domain 1
bndMesh1 = [];
for i=1:mesh1.noElemsV
    bndMesh1 = [bndMesh1; mesh1.noElemsU*i ];
end

% boundary edges in V-direction
eps = 1e-14;
bndNodes  = find(abs(mesh1.controlPts(:,1)-L1) < eps);
bndEdge1  = zeros(mesh1.noElemsV,mesh1.q+1);

for i=1:mesh1.noElemsV
    bndEdge1(i,:) = bndNodes(i:i+mesh1.q);
end

%% Gauss points on coupling interface
ngp = 12;
[W1,Q1]=quadrature( ngp, 'GAUSS', 1 ); % two point quadrature

GP1 = [];
GP2 = [];

for e=1:mesh1.noElemsV
    sctrEdge = bndEdge1(e,:);
    sctr1    = mesh1.globElems(bndMesh1(e),:);
    pts1     = mesh1.controlPts(sctr1,:);
    ptsE     = mesh1.controlPts(sctrEdge,:);
    xiE1     = mesh1.rangeU(end,:); % [xi_i,xi_i+1]
    etE1     = mesh1.rangeV(e,:);   % [et_j,eta_j+1]
    xiE2     = mesh2.range(1,:);   % [et_j,eta_j+1]
    for q=1:size(W1)
        pt=Q1(q,:);
        wt=W1(q);
        Xi      = 0.5 * ( ( etE1(2) - etE1(1) ) * pt + etE1(2) + etE1(1));
        J2      = 0.5 * ( etE1(2) - etE1(1) );
        [N, dNdxi] = NURBS1DBasisDers(Xi,mesh1.q,mesh1.vKnot,mesh1.weights);
        J0    = dNdxi*ptsE; % also the tangent to Gamma_*^e
        detJ0 = norm(J0);
        J0    = J0/detJ0;
        
        x  = N*ptsE;
        X1 = global2LocalMapNURBS2D(x,pts1,xiE1,etE1,mesh1.p,mesh1.q, ...
                                    mesh1.uKnot, mesh1.vKnot, mesh1.weights);
        X2  = 0.5 * ( ( xiE2(2) - xiE2(1) ) * (-1) + xiE2(2) + xiE2(1));
        GP1 = [GP1;X1 wt*detJ0*J2 -J0(2) J0(1)];
        GP2 = [GP2;X2 x(2)];
    end
end

%% DEFINE BOUNDARIES

% Dirichlet nodes and other node groups
fixedNode = find(abs(mesh1.controlPts(:,1)) < eps);
midNode1  = find(abs(mesh1.controlPts(:,2)) < eps);

uFixed=zeros(1,length(fixedNode))';  % a vector of u_x for the nodes
vFixed=zeros(1,length(fixedNode))';


%% least square for boundary conditions

dispNodes   = fixedNode;
noDispNodes = length(dispNodes);


bndElement     = zeros(mesh1.noElemsV,mesh1.q+1); % assume p=q!!!
leftEdgeMesh   = zeros(mesh1.noElemsV,mesh1.q+1);

for i=1:mesh1.noElemsV
    leftEdgeMesh(i,:)  = fixedNode (i:i+mesh1.q);
end

% build boundary mesh for assembly matrix A in Aq=b

for ie=1:mesh1.noElemsV
    bndElement(ie,:) = [ie:ie+mesh1.q];
end

% essential boundary conditions, exact displacement is used
% on the boundary edges

disp([num2str(toc),'  LEAST SQUARE for DIRICHLET BCs'])

A  = zeros(noDispNodes,noDispNodes);
bx = zeros(noDispNodes,1);
by = zeros(noDispNodes,1);

noxC   = 16;

% Loop over boundary edges...
% left edge

for ie=1:mesh1.noElemsV
    sctr   = leftEdgeMesh(ie,:);
    pts    = mesh1.controlPts(sctr,:);
    sctrA  = bndElement(ie,:);
    xiE    = mesh1.rangeV(ie,:);
    xiArr  = linspace(xiE(1),xiE(2),noxC);
    
    for ic=1:noxC
        xi = xiArr(ic);
        [N dNdxi] = NURBS1DBasisDers(xi,mesh1.q,mesh1.vKnot,mesh1.weights);
        
        A(sctrA,sctrA) = A(sctrA,sctrA) + N'*N;
        
        pt= N    *pts;
        x = pt(1); y = pt(2);
        
        % exact displacements
        
        ux = 1000*y/(6*E0*I0)*((6*L-3*x)*x+(2+nu0)*(y^2-c^2));
        uy = -1000/(6*E0*I0)*(3*nu0*y^2*(L-x)+(4+5*nu0)*t^2*x/4+(3*L-x)*x^2);
        
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

uFixed     = qx;
vFixed     = qy;

% for i=1:length(fixedNode)
%     inode = fixedNode(i);
%     pts   = mesh1.controlPts(inode,:);
%     ux    = 1000*pts(2)/(6*E0*I0)*(2+nu0)*(pts(2)^2-c^2);
%     uy    = -1000/(2*E0*I0)*nu0*pts(2)^2*L;
%     uFixed(i)=ux;
%     vFixed(i)=uy;
% end

%%
numnodes = size(mesh1.controlPts,1) + size(mesh2.controlPts,1);

fiNodes=[fixedNode;numnodes];
frNodes=setdiff(1:numnodes,fiNodes);

activeDofs=[];
for i=1:length(frNodes)
    in = frNodes(i);
    activeDofs=[activeDofs; 2*in-1;2*in];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE SYSTEM DATA STRUCTURES

disp([num2str(toc),' INITIALIZING DATA STRUCTURES'])

f=zeros(2*numnodes,1);          % external load vector
K=zeros(2*numnodes,2*numnodes); % stiffness matrix

% ******************************************************************************
% ***                          P R O C E S S I N G                           ***
% ******************************************************************************
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%% COMPUTE STIFFNESS MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp([num2str(toc),'   COMPUTING STIFFNESS MATRIX'])

%% for domain 1 -----------------------------------------------------------

ngpv = mesh1.p+1;

[W,Q]=quadrature( ngpv, 'GAUSS', 2 ); % 2x2 Gaussian quadrature

for e=1:mesh1.elemCount                    % loop over elements
    idu    = mesh1.index(e,1);
    idv    = mesh1.index(e,2);
    xiE    = mesh1.rangeU(idu,:); % [xi_i,xi_i+1]
    etaE   = mesh1.rangeV(idv,:); % [eta_j,eta_j+1]
    
    sctr   = mesh1.globElems(e,:);         %  element scatter vector
    nn     = length(sctr);
    
    sctrB(1,1:2:2*nn) = 2*sctr-1;
    sctrB(1,2:2:2*nn) = 2*sctr  ;
    
    pts    = mesh1.controlPts(sctr,:);
      
    for gp=1:size(W,1)                     % loop over Gauss points
        pt      = Q(gp,:);
        wt      = W(gp);
        Xi      = parent2ParametricSpace(xiE, pt(1));
        Eta     = parent2ParametricSpace(etaE,pt(2));
        J2      = jacobianPaPaMapping(xiE,etaE);
        
        % compute derivatives of shape functions
        [R, dRdxi, dRdeta] = NURBS2DBasisDers([Xi;Eta],mesh1.p,mesh1.q,mesh1.uKnot,mesh1.vKnot,mesh1.weights');
        
        % compute the jacobian of physical and parameter domain mapping
        % then the derivative w.r.t spatial physical coordinates
        
        jacob      = pts' * [dRdxi' dRdeta'];
        J1         = det(jacob);
        dRdx       = [dRdxi' dRdeta'] * inv(jacob);
        J          = J1 * J2;
        B           = getBmatrix2D(dRdx');
        K(sctrB,sctrB) = K(sctrB,sctrB) + B' * C * B * J * wt;
    end
end

%% for domain 2 (beam) ----------------------------------------------------
clear sctrB;
noGPs=solid2.order;
[W,Q]=quadrature(  noGPs, 'GAUSS', 1 ); % 2 point quadrature

for e=1:mesh2.elemCount
   xiE   = mesh2.range(e,:); % [xi_i,xi_i+1]
   conn  = mesh2.element(e,:);
   noFns = length(conn);
   conng = conn + size(mesh1.controlPts,1);
   sctrB(1:2:2*noFns) = 2*conng-1;
   sctrB(2:2:2*noFns) = 2*conng-0;
   
   pts   = mesh2.controlPts(conn,1:2);
      
    for gp=1:size(W,1)                 % loop over Gauss points                  
      pt      = Q(gp,:);                          
      wt      = W(gp);                 
      % coord in parameter space  
      Xi      = 0.5 * (( xiE(2) - xiE(1) ) * pt + xiE(2) + xiE(1));
      J2      = 0.5 * ( xiE(2) - xiE(1) ); 
      
      % shape functions and first derivatives             
      [R, dRdxi] = NURBS1DBasisDers(Xi,solid2.order-1,solid2.knots,solid2.coefs(4,:));
      
      % compute the jacobian of physical and parameter domain mapping
      % then the derivative w.r.t spatial physical coordinates
      
      dxdxi  = dRdxi *pts;      
      J1     = norm (dxdxi);
      dRdx   = (1/J1)*dRdxi;
                
      Bb(2:2:noFns*2)   = dRdx;
      Bs(1:2:noFns*2)   = dRdx;
      Bs(2:2:noFns*2)   = -R;      
      % compute elementary stiffness matrix and
      % assemble it to the global matrix
      K(sctrB,sctrB) = K(sctrB,sctrB) + ...
          (k * Bb' * Bb + shearFactor*G*b*t*Bs'*Bs)* J1 * J2 * wt;
    end
end
clear sctrB;

%% determination of alpha

Ktilde = K;
H=zeros(2*numnodes,2*numnodes); % stiffness matrix

Cbeam  = [E0 0;0 shearFactor*G];
Csolid = C([1 3],:);

for i=1:mesh1.noElemsV                     % start of element loop
    e1     = bndMesh1(i);
    sctr1  = mesh1.globElems(e1,:);
    sctr2  = mesh2.element(1,:);
    sctr2n = sctr2 + size(mesh1.controlPts,1);
    nn1    = length(sctr1);
    nn2    = length(sctr2);
    sctrB1(1:2:2*nn1) = 2*sctr1-1;
    sctrB1(2:2:2*nn1) = 2*sctr1-0;
    sctrB2(1:2:2*nn2) = 2*sctr2n-1;
    sctrB2(2:2:2*nn2) = 2*sctr2n-0;
    
    pts1 = mesh1.controlPts(sctr1,:);
    pts2 = mesh2.controlPts(sctr2,:);
        
    Kd11 = zeros(nn1*2,nn1*2);
    Kd12 = zeros(nn1*2,nn2*2);
    Kd21 = zeros(nn2*2,nn1*2);
    Kd22 = zeros(nn2*2,nn2*2);
    
    for q=ngp*(i-1)+1:ngp*(i-1)+ngp        
        pt1=GP1(q,1:2);
        wt1=GP1(q,3);
        pt2=GP2(q,1);
        y  =GP2(q,2);
        normal=GP1(q,4:5);
        n = [normal(1) normal(2);0 normal(1)];
        n=-n;
        [N1 dN1dxi dN1deta] = NURBS2DBasisDers(pt1,mesh1.p,mesh1.q,mesh1.uKnot,mesh1.vKnot,mesh1.weights');
        [N2 dN2dxi]         = NURBS1DBasisDers(pt2,solid2.order-1,solid2.knots,solid2.coefs(4,:));
        
        J1 = pts1' * [dN1dxi' dN1deta'];
        J2 = dN2dxi*pts2;
        detJ0=norm(J2);
        invJacob1   = inv(J1);
        dN1dx      = [dN1dxi' dN1deta'] * invJacob1;
        
        dN2dx   = (1/detJ0)*dN2dxi;
        
        B1(1,1:2:2*nn1)  = dN1dx(:,1)';
        B1(2,2:2:2*nn1)  = dN1dx(:,2)';
        B1(3,1:2:2*nn1)  = dN1dx(:,2)';
        B1(3,2:2:2*nn1)  = dN1dx(:,1)';
        
        B2(1,2:2:2*nn2)  = -y*dN2dx(:,1)';
        B2(2,1:2:2*nn2)  = dN2dx(:,1)';
        B2(2,2:2:2*nn2)  = -N2;
        
        Nm1(1,1:2:2*nn1)  = N1';
        Nm1(2,2:2:2*nn1)  = N1';
        
        Nm2(1,2:2:2*nn2)  = -y*N2';
        Nm2(2,1:2:2*nn2)  = N2';
        
        Kd11 = Kd11 + B1'*Csolid'*n'*n*Csolid*B1 *wt1;
        Kd12 = Kd12 + B1'*Csolid'*n'*n*Cbeam*B2  *wt1;
        Kd21 = Kd21 + B2'*Cbeam' *n'*n*Csolid*B1 *wt1;
        Kd22 = Kd22 + B2'*Cbeam' *n'*n*Cbeam*B2  *wt1;  
        
    end  % of quadrature loop
    
    H(sctrB1,sctrB1)  = H(sctrB1,sctrB1)  + Kd11;
    H(sctrB1,sctrB2)  = H(sctrB1,sctrB2)  + Kd12;
    H(sctrB2,sctrB1)  = H(sctrB2,sctrB1)  + Kd21;
    H(sctrB2,sctrB2)  = H(sctrB2,sctrB2)  + Kd22;
end

lambda = eigs(inv(Ktilde(activeDofs,activeDofs))*H(activeDofs,activeDofs),5,'lm');

alpha=max(lambda)/2;

%alpha=3e8;

%% interface integrals ----------------------------------------------------

Cbeam  = [E0 0;0 G];
Csolid = C([1 3],:);
beta=1/c^2;
beta=0;

for i=1:mesh1.noElemsV                     % start of element loop
    e1     = bndMesh1(i);
    sctr1  = mesh1.globElems(e1,:);
    sctr2  = mesh2.element(1,:);
    sctr2n = sctr2 + size(mesh1.controlPts,1);
    nn1    = length(sctr1);
    nn2    = length(sctr2);
    sctrB1(1:2:2*nn1) = 2*sctr1-1;
    sctrB1(2:2:2*nn1) = 2*sctr1-0;
    sctrB2(1:2:2*nn2) = 2*sctr2n-1;
    sctrB2(2:2:2*nn2) = 2*sctr2n-0;
    
    pts1 = mesh1.controlPts(sctr1,:);
    pts2 = mesh2.controlPts(sctr2,:);
    
    Kp11 = zeros(nn1*2,nn1*2);
    Kp12 = zeros(nn1*2,nn2*2);
    Kp22 = zeros(nn2*2,nn2*2);
    
    Kd11 = zeros(nn1*2,nn1*2);
    Kd12 = zeros(nn1*2,nn2*2);
    Kd21 = zeros(nn2*2,nn1*2);
    Kd22 = zeros(nn2*2,nn2*2);
    
    for q=ngp*(i-1)+1:ngp*(i-1)+ngp        
        pt1 = GP1(q,1:2);
        wt1 = GP1(q,3);
        pt2 = GP2(q,1);
        y   = GP2(q,2);
        normal=GP1(q,4:5);
        n = [normal(1) normal(2);0 normal(1)];
        n=-n;
        [N1, dN1dxi, dN1deta] = NURBS2DBasisDers(pt1,mesh1.p,mesh1.q,mesh1.uKnot,mesh1.vKnot,mesh1.weights');
        [N2, dN2dxi, dN2dxi2] = NURBS1DBasis2ndDers(pt2,solid2.order-1,solid2.knots,solid2.coefs(4,:));
                    
        J2      = dN2dxi *pts2;
        dx2dxi2 = dN2dxi2*pts2;
        detJ0     = norm(J2);
        dN2dx     = (1/detJ0)*dN2dxi;                
        dN2dx2 = 1/detJ0^2 * (dN2dxi - norm(dx2dxi2)*dN2dx);
      
        J1        = pts1' * [dN1dxi' dN1deta'];                        
        dN1dx     = [dN1dxi' dN1deta'] * inv(J1);        
        
        
        B1(1,1:2:2*nn1)  = dN1dx(:,1)';
        B1(2,2:2:2*nn1)  = dN1dx(:,2)';
        B1(3,1:2:2*nn1)  = dN1dx(:,2)';
        B1(3,2:2:2*nn1)  = dN1dx(:,1)';
                
        Nm1(1,1:2:2*nn1)  = N1;
        Nm1(2,2:2:2*nn1)  = N1;
        
        B2(1,1:2:2*nn2)  = (-beta/3*y^3) *dN2dx2;
        B2(1,2:2:2*nn2)  = (-y+beta/3*y^3)*dN2dx;
        B2(2,1:2:2*nn2)  = (1-beta*y^2)  *dN2dx;
        B2(2,2:2:2*nn2)  = (1-beta*y^2)  *N2;
                       
        Nm2(1,1:2:2*nn2)  =  (-beta/3*y^3) *dN2dx;
        Nm2(1,2:2:2*nn2)  =  (-y+beta/3*y^3)*N2;
        Nm2(2,1:2:2*nn2)  = (1+y^2/c^2)*N2;
        
        %Cbeam(2,2) = 5*E0/8/(1+nu0);
        
        % COMPUTE ELEMENT STIFFNESS AT QUADRATURE POINT
        
        Kp11 = Kp11 + alpha*(Nm1'*Nm1)*wt1;
        Kp12 = Kp12 + alpha*(Nm1'*Nm2)*wt1;
        Kp22 = Kp22 + alpha*(Nm2'*Nm2)*wt1;
        
        Kd11 = Kd11 + gamma1 * Nm1'* n * Csolid * B1 *wt1;
        Kd12 = Kd12 + gamma2 * Nm1'* n * Cbeam  * B2 *wt1;
        Kd21 = Kd21 + gamma1 * Nm2'* n * Csolid * B1 *wt1;
        Kd22 = Kd22 + gamma2 * Nm2'* n * Cbeam  * B2 *wt1;
        
    end  % of quadrature loop
    
    K(sctrB1,sctrB1)  = K(sctrB1,sctrB1)  - Kd11 - Kd11' + Kp11;
    K(sctrB1,sctrB2)  = K(sctrB1,sctrB2)  - Kd12 + Kd21' - Kp12;
    K(sctrB2,sctrB1)  = K(sctrB2,sctrB1)  + Kd21 - Kd12' - Kp12';
    K(sctrB2,sctrB2)  = K(sctrB2,sctrB2)  + Kd22 + Kd22' + Kp22;
end


f(2*numnodes-1) = P;

%%%%%%%%%%%%%%%%%%% END OF STIFFNESS MATRIX COMPUTATION %%%%%%%%%%%%%%%%%%%


% APPLY ESSENTIAL BOUNDARY CONDITIONS
disp([num2str(toc),'   APPLYING BOUNDARY CONDITIONS'])
bcwt=mean(diag(K)); % a measure of the average size of an element in K
% used to keep the conditioning of the K matrix
udofs=2*fixedNode-1;           % global indecies of the fixed x displacements
vdofs=2*fixedNode;   % global indecies of the fixed y displacements
f=f-K(:,udofs)*uFixed;  % modify the force vector
f=f-K(:,vdofs)*vFixed;
K(udofs,:)=0;
K(vdofs,:)=0;
K(:,udofs)=0;
K(:,vdofs)=0;
K(udofs,udofs)=bcwt*speye(length(udofs)); % put ones*bcwt on the diagonal
K(vdofs,vdofs)=bcwt*speye(length(vdofs));
f(udofs)=bcwt*speye(length(udofs))*uFixed;
f(vdofs)=bcwt*speye(length(vdofs))*vFixed;


% SOLVE SYSTEM
disp([num2str(toc),'   SOLVING SYSTEM'])
U=K\f;

%******************************************************************************
%*** POST - PROCESSING ***
%***************************************************

disp([num2str(toc),'   POST-PROCESSING'])

numnode1 = size(mesh1.controlPts,1);
numnode2 = size(controlPts,1);

Ux = U(1:2:2*numnodes);
Uy = U(2:2:2*numnodes);

Ux1 = Ux(1:numnode1);
Ux2 = Ux(numnode1+1:end);

Uy1 = Uy(1:numnode1);
Uy2 = Uy(numnode1+1:end);

%%

Uxm1  = U(2*midNode1-1);
Uym1  = U(2*midNode1);


xx = mesh1.controlPts(midNode1,1);
u  = Uym1;

xx = [xx;mesh2.controlPts(:,1)];
u  = [u;Ux2];

% exact solution u =
y= 0;
x      = linspace(0,L,200);
D=t;
uExact = -1000/6/E0/I0*(3*nu0*y*y*(L-x)+(4+5*nu0)*D*D*x/4+(3*L-x).*x.^2);

colordef white
figure,set (gcf,'Color','w')
set(gca,'FontSize',14)
hold on
plot(x,uExact,'k-','LineWidth',1.4);
plot(xx,u,'o','MarkerEdgeColor','k',...
    'MarkerFaceColor','g',...
    'MarkerSize',6.5);
h=legend('exact','coupling');
xlabel('x')
ylabel('w')
grid on

%% plot stresses at midline

aa=floor(mesh1.noElemsV/2);
aa=mesh1.noElemsV-1;
elems1 = mesh1.noElemsU*aa+1:mesh1.noElemsU*(aa+1);


sigma      = zeros(length(elems1),3);
sigmaRef   = zeros(length(elems1),3);
xcoord     = zeros(length(elems1),2);

for e=1:length(elems1)
    ie = elems1(e);
    idu    = mesh1.index(ie,1);
    idv    = mesh1.index(ie,2);
    xiE    = mesh1.rangeU(idu,:); % [xi_i,xi_i+1]
    etaE   = mesh1.rangeV(idv,:); % [eta_j,eta_j+1]
    
    sctr   = mesh1.globElems(ie,:);         %  element scatter vector
    nn     = length(sctr);
    
    sctrB(1,1:2:2*nn) = 2*sctr-1;
    sctrB(1,2:2:2*nn) = 2*sctr  ;
    
    pts    = mesh1.controlPts(sctr,:);
    
    % loop over Gauss points
    
    %for gp=1:size(W,1)
        pt      = [0 0];
        Xi      = parent2ParametricSpace(xiE,pt(1));
        Eta     = parent2ParametricSpace(etaE,pt(2));
        J2      = jacobianPaPaMapping(xiE,etaE);
        
        % compute derivatives of shape functions
        [R, dRdxi, dRdeta] = NURBS2DBasisDers([Xi;Eta],mesh1.p,mesh1.q,mesh1.uKnot,mesh1.vKnot,mesh1.weights');
        
        % compute the jacobian of physical and parameter domain mapping
        % then the derivative w.r.t spatial physical coordinates
        
        jacob      =  pts' * [dRdxi' dRdeta'];        
        invJacob   = inv(jacob);
        dRdx       = [dRdxi' dRdeta'] * invJacob;  
        yPt=R*pts;
        B           = getBmatrix2D(dRdx');
        % COMPUTE ELEMENT STRAIN AND STRESS AT STRESS POINT
        strain=B*U(sctrB);
        stress=C*strain;
        sigma(e,1)    = stress(1);
        sigma(e,2)    = stress(3);
        sigma(e,3)    = stress(2);
        sigmaRef(e,1) = 1000/I0*(L-yPt(1))*yPt(2);
        sigmaRef(e,2) = -1000/2/I0*(c^2-yPt(2)^2);
        sigmaRef(e,3) = 0;
        xcoord(e,:)     = yPt;
end   % of element loop
    
colordef white
figure,set (gcf,'Color','w')
set(gca,'FontSize',14)
hold on
% plot(xcoord(:,1),sigmaRef(:,1),'k-','LineWidth',1.4);
% plot(xcoord(:,1),sigma(:,1),'o','MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',6.5);
%plot(xcoord(:,1),sigmaRef(:,2),'k--','LineWidth',1.4);
%plot(xcoord(:,1),sigma(:,2),'s','MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',6.5);
plot(xcoord(:,1),sigmaRef(:,3),'k-','LineWidth',1.4);
plot(xcoord(:,1),sigma(:,3),'o','MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',6.5);
h=legend('sigmaxx-exact','sigmaxx-coupling','sigmaxy-exact','sigmaxy-coupling');
xlabel('x')
ylabel('stresses along y=0.3')
grid on
    
%%
    
aa=mesh1.noElemsU-0;
elems1 = aa:mesh1.noElemsU:aa+mesh1.noElemsU*(mesh1.noElemsV-1);

sigma      = zeros(length(elems1),2);
sigmaRef   = zeros(length(elems1),2);
xcoord     = zeros(length(elems1),2);

for e=1:length(elems1)
    ie = elems1(e);
    idu    = mesh1.index(ie,1);
    idv    = mesh1.index(ie,2);
    xiE    = mesh1.rangeU(idu,:); % [xi_i,xi_i+1]
    etaE   = mesh1.rangeV(idv,:); % [eta_j,eta_j+1]
    
    sctr   = mesh1.globElems(ie,:);         %  element scatter vector
    nn     = length(sctr);
    
    sctrB(1,1:2:2*nn) = 2*sctr-1;
    sctrB(1,2:2:2*nn) = 2*sctr  ;
    
    pts    = mesh1.controlPts(sctr,:);
    
    % loop over Gauss points
    
    %for gp=1:size(W,1)
    pt      = [-0 0];
    Xi      = parent2ParametricSpace(xiE,pt(1));
    Eta     = parent2ParametricSpace(etaE,pt(2));
    J2      = jacobianPaPaMapping(xiE,etaE);
    
    % compute derivatives of shape functions
    [R dRdxi dRdeta] = NURBS2DBasisDers([Xi; Eta],mesh1.p,mesh1.q,mesh1.uKnot,mesh1.vKnot,mesh1.weights');
    
    % compute the jacobian of physical and parameter domain mapping
    % then the derivative w.r.t spatial physical coordinates
    
    jacob      =  pts' * [dRdxi' dRdeta'];
    invJacob   = inv(jacob);
    dRdx       = [dRdxi' dRdeta'] * invJacob;
    yPt=R*pts;
    B           = getBmatrix2D(dRdx');
    % COMPUTE ELEMENT STRAIN AND STRESS AT STRESS POINT
    strain=B*U(sctrB);
    stress=C*strain;
    sigma(e,1)    = stress(1);
    sigma(e,2)    = stress(3);
    sigmaRef(e,1) = 1000/I0*(L-yPt(1))*yPt(2);
    sigmaRef(e,2) = -1000/2/I0*(c^2-yPt(2)^2);
    xcoord(e,:)     = yPt;
end   % of element loop


colordef white
figure,set (gcf,'Color','w')
set(gca,'FontSize',14)
hold on
% plot(xcoord(:,2),sigmaRef(:,1),'k-','LineWidth',1.4);
% plot(xcoord(:,2),sigma(:,1),'o','MarkerEdgeColor','k',...
%     'MarkerFaceColor','g',...
%     'MarkerSize',6.5);
plot(xcoord(:,2),sigmaRef(:,2),'k-','LineWidth',1.4);
plot(xcoord(:,2),sigma(:,2),'s','MarkerEdgeColor','k',...
    'MarkerFaceColor','g',...
    'MarkerSize',6.5);
h=legend('sigmaxy-exact','sigmaxy-coupling');
xlabel('y')
ylabel('sigma_{xy} at x=L/2')
grid on

%% plot stresses on coupling interface

%aa=mesh1.noElemsU;
%elems1 = aa:mesh1.noElemsU:aa+mesh1.noElemsU*(mesh1.noElemsV-1);

% sigma      = zeros(length(elems1),2);
% sigmaRef   = zeros(length(elems1),2);
% xcoord     = zeros(length(elems1),2);
% 
% e = 1;
% for i=1:mesh1.noElemsV                     % start of element loop
%     e1     = bndMesh1(i);
%     sctr1  = mesh1.globElems(e1,:);        
%     nn1    = length(sctr1);    
%     sctrB1(1:2:2*nn1) = 2*sctr1-1;
%     sctrB1(2:2:2*nn1) = 2*sctr1-0;
%     
%     pts1 = mesh1.controlPts(sctr1,:);
% 
%     for q=ngp*(i-1)+1:ngp*(i-1)+ngp        
%         pt1 = GP1(q,1:2);
%         wt1 = GP1(q,3);
%         [N1, dN1dxi, dN1deta] = NURBS2DBasisDers(pt1,mesh1.p,mesh1.q,mesh1.uKnot,mesh1.vKnot,mesh1.weights');
%                 
%         J1        = pts1' * [dN1dxi' dN1deta'];                        
%         dN1dx     = [dN1dxi' dN1deta'] * inv(J1);                
%         
%         B1(1,1:2:2*nn1)  = dN1dx(:,1)';
%         B1(2,2:2:2*nn1)  = dN1dx(:,2)';
%         B1(3,1:2:2*nn1)  = dN1dx(:,2)';
%         B1(3,2:2:2*nn1)  = dN1dx(:,1)';   
%         
%         yPt=N1*pts1;
%         
%         % COMPUTE ELEMENT STRAIN AND STRESS AT STRESS POINT
%         strain=B1*U(sctrB1);
%         stress=C*strain;
%         sigma(e,1)    = stress(1);
%         sigma(e,2)    = stress(3);
%         sigmaRef(e,1) = 1000/I0*(L-yPt(1))*yPt(2);
%         sigmaRef(e,2) = -1000/2/I0*(c^2-yPt(2)^2);
%         xcoord(e,:)     = yPt;
%     
%         e = e + 1;
%     end  % of quadrature loop
% end
% 

aa=mesh1.noElemsU-0;
elems1 = aa:mesh1.noElemsU:aa+mesh1.noElemsU*(mesh1.noElemsV-1);

disp       = zeros(length(elems1),2);
sigma      = zeros(length(elems1),2);
dispRef    = zeros(length(elems1),2);
sigmaRef   = zeros(length(elems1),2);
xcoord     = zeros(length(elems1),2);
err1=0;
err2=0;
err3=0;
for e=1:length(elems1)
    ie = elems1(e);
    idu    = mesh1.index(ie,1);
    idv    = mesh1.index(ie,2);
    xiE    = mesh1.rangeU(idu,:); % [xi_i,xi_i+1]
    etaE   = mesh1.rangeV(idv,:); % [eta_j,eta_j+1]
    
    sctr   = mesh1.globElems(ie,:);         %  element scatter vector
    nn     = length(sctr);
    
    sctrB(1,1:2:2*nn) = 2*sctr-1;
    sctrB(1,2:2:2*nn) = 2*sctr  ;
    sctry = 2*sctr;
    
    pts    = mesh1.controlPts(sctr,:);
    
    % loop over Gauss points
    
    %for gp=1:size(W,1)
    pt      = [1 0];
    Xi      = parent2ParametricSpace(xiE,pt(1));
    Eta     = parent2ParametricSpace(etaE,pt(2));
    J2      = jacobianPaPaMapping(xiE,etaE);
    
    % compute derivatives of shape functions
    [R, dRdxi, dRdeta] = NURBS2DBasisDers([Xi;Eta],mesh1.p,mesh1.q,mesh1.uKnot,mesh1.vKnot,mesh1.weights');
    
    % compute the jacobian of physical and parameter domain mapping
    % then the derivative w.r.t spatial physical coordinates
    
    jacob      =  pts' * [dRdxi' dRdeta'];
    invJacob   = inv(jacob);
    dRdx       = [dRdxi' dRdeta'] * invJacob;
    yPt=R*pts;
    B           = getBmatrix2D(dRdx');
    % COMPUTE ELEMENT STRAIN AND STRESS AT STRESS POINT
    strain=B*U(sctrB);
    stress=C*strain;
    disp(e,2)     = R*U(sctry);
    sigma(e,1)    = stress(1);
    sigma(e,2)    = stress(3);
    sigmaRef(e,1) = 1000/I0*(L-yPt(1))*yPt(2);
    sigmaRef(e,2) = -1000/2/I0*(c^2-yPt(2)^2);
    dispRef(e,1)  = -1000/6/E0/I0*(3*nu0*yPt(2)^2*(L-yPt(1))+(4+5*nu0)*D*D*yPt(1)/4+(3*L-yPt(1))*yPt(1)^2);
    xcoord(e,:)   = yPt;
    der1 = (sigmaRef(e,1)-sigma(e,1))/sigmaRef(e,1);
    der2 = (sigmaRef(e,2)-sigma(e,2))/sigmaRef(e,2);
    der3 = (dispRef(e,1) -disp)/dispRef(e,1);
    err1 = err1 + der1^2;
    err2 = err2 + der2^2;
%    err3 = err3 + der3^2;
end   % of element loop

err3 = sqrt(err3) % error in disp (vertical)
err1 = sqrt(err1) % error in sigma_xx
err2 = sqrt(err2) % error in sigma_xy


% 
% % 
colordef white
figure,set (gcf,'Color','w')
set(gca,'FontSize',14)
hold on
plot(xcoord(:,2),disp(:,2),'k-','LineWidth',1.4);
%plot(xcoord(:,2),sigma(:,2),'o','MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',6.5);
% plot(xcoord(:,2),sigmaRef(:,2),'k-','LineWidth',1.4);
% plot(xcoord(:,2),sigma(:,2),'s','MarkerEdgeColor','k',...
%     'MarkerFaceColor','g',...
%     'MarkerSize',6.5);
h=legend('sigmaxy-exact','sigmaxy-coupling');
xlabel('y')
ylabel('sigma_{xy} at x=L/2')
grid on

%% plot contour of stress field

material.stiffMat=C;
displacement = [U(1:2:end) U(2:2:end)];

matMap1 = ones(mesh1.elemCount,1);

vMesh1=buildVisualizationMesh2D(solid1);

meshes{1}  = mesh1;
vmeshes{1} = vMesh1;

meshData.mesh    = meshes;
meshData.vmesh   = vmeshes;
meshData.matMap{1}=matMap1;
materials{1}      = material;
meshData.voids   = [];

vtuFileName = 'iga-nitsche-2D-beam';
for ip=1:1
    vtuFile = strcat(vtuFileName,'-mesh',num2str(ip));
    %figure; hold on;
    ok      = plotStress2DForPatch(meshData,ip,vtuFile,displacement,materials);
end
% 
% Cinv           = inv(C);
% energy_norm    = 0;
% disp_norm      = 0;
% 
% 
% for ip=1:length(meshes)
%     mesh = meshes{ip};
%     [W,Q]= quadrature( mesh.p+1, 'GAUSS', 2 ); % 2x2 Gaussian quadrature
%     
%     for e=1:mesh.elemCount
%         idu    = mesh.index(e,1);
%         idv    = mesh.index(e,2);
%         xiE    = mesh.rangeU(idu,:); % [xi_i,xi_i+1]
%         etaE   = mesh.rangeV(idv,:); % [eta_j,eta_j+1]
%         
%         sctrg  = mesh.globElems(e,:);         %  element scatter vector
%         sctrl  = mesh.locElems (e,:);         %  element scatter vector
%         sctrx  = 2*sctrg-1;
%         sctry  = 2*sctrg  ;
%         nn     = length(sctrg);
%         sctrB  = zeros(1,2*nn);
%         sctrB(1,1:2:2*nn) = 2*sctrg-1;
%         sctrB(1,2:2:2*nn) = 2*sctrg  ;
%         
%         pts    = mesh.controlPts(sctrl,:);
%         
%         % loop over Gauss points
%         
%         for gp=1:size(W,1)
%             pt      = Q(gp,:);
%             wt      = W(gp);
%             
%             % compute coords in parameter space
%             Xi      = parent2ParametricSpace(xiE, pt(1));
%             Eta     = parent2ParametricSpace(etaE,pt(2));
%             J2      = jacobianPaPaMapping(xiE,etaE);
%             % compute derivatives of shape functions
%             [R dRdxi dRdeta] = NURBS2DBasisDers([Xi;Eta],mesh.p,mesh.q,mesh.uKnot,mesh.vKnot,mesh.weights');
%             
%             % compute the jacobian of physical and parameter domain mapping
%             % then the derivative w.r.t spatial physical coordinates
%             
%             jacob      = pts' * [dRdxi' dRdeta'];
%             J1         = det(jacob);
%             invJacob   = inv(jacob);
%             dRdx       = [dRdxi' dRdeta'] * invJacob;
%             J          = J1 * J2;
%             B          = getBmatrix2D(dRdx');
%             % Numerical strain and stress and displacements            
%             strain     = B*U(sctrB);
%             num_stress = C*strain;
%             num_disp   = R*[U(sctrx) U(sctry)];
%             strPoint   = R * pts;
%             % Exact stresses and displacements
%             
%             [exact_stress, exact_disp]= exact_solution_timoshenkobeam(...
%                 strPoint,E0,nu0,I0,L,2*c,1000);
%             
%             errorInNorm = num_stress - exact_stress';
%             errorInDisp = num_disp   - exact_disp;
%             
%             fac            = wt*J1*J2;
%             energy_norm    = energy_norm    + (Cinv*errorInNorm)'*errorInNorm*fac;            
%             disp_norm      = disp_norm      + errorInDisp*errorInDisp'*fac;
%         end
%     end
% end
% 
% disp_norm   = sqrt(disp_norm)
% energy_norm = sqrt(0.5*energy_norm)
% 2*numnodes
Cbeam  = [E0 0;0 shearFactor*G];
sigmaSolid  = zeros(length(elems1),2);
sigmaBeam   = zeros(length(elems1),2);
sigmaRef    = zeros(length(elems1),2);
xcoord      = zeros(length(elems1),2);
dispRef     = zeros(length(elems1),2);
dispSolid   = zeros(length(elems1),2);
dispBeam    = zeros(length(elems1),2);

e = 1;
for i=1:mesh1.noElemsV                     % start of element loop
    e1     = bndMesh1(i);
    sctr1  = mesh1.globElems(e1,:);
    sctr2  = mesh2.element(1,:);
    sctr2n = sctr2 + size(mesh1.controlPts,1);
    nn1    = length(sctr1);
    nn2    = length(sctr2);
    sctrB1(1:2:2*nn1) = 2*sctr1-1;
    sctrB1(2:2:2*nn1) = 2*sctr1-0;
    sctrB2(1:2:2*nn2) = 2*sctr2n-1;
    sctrB2(2:2:2*nn2) = 2*sctr2n-0;
    
    pts1 = mesh1.controlPts(sctr1,:);
    pts2 = mesh2.controlPts(sctr2,:);

    for q=ngp*(i-1)+1:ngp*(i-1)+ngp        
        pt1 = GP1(q,1:2);
        wt1 = GP1(q,3);
        pt2 = GP2(q,1);
        y   = GP2(q,2);
        normal=GP1(q,4:5);
        n = [normal(1) normal(2);0 normal(1)];
        n=-n;
        [N1, dN1dxi, dN1deta] = NURBS2DBasisDers(pt1,mesh1.p,mesh1.q,mesh1.uKnot,mesh1.vKnot,mesh1.weights');
        [N2, dN2dxi, dN2dxi2] = NURBS1DBasis2ndDers(pt2,solid2.order-1,solid2.knots,solid2.coefs(4,:));
                    
        J2      = dN2dxi *pts2;
        dx2dxi2 = dN2dxi2*pts2;
        detJ0     = norm(J2);
        dN2dx     = (1/detJ0)*dN2dxi;                
        dN2dx2 = 1/detJ0^2 * (dN2dxi - norm(dx2dxi2)*dN2dx);
      
        J1        = pts1' * [dN1dxi' dN1deta'];                        
        dN1dx     = [dN1dxi' dN1deta'] * inv(J1);        
        
        
        B1(1,1:2:2*nn1)  = dN1dx(:,1)';
        B1(2,2:2:2*nn1)  = dN1dx(:,2)';
        B1(3,1:2:2*nn1)  = dN1dx(:,2)';
        B1(3,2:2:2*nn1)  = dN1dx(:,1)';
                
        Nm1(1,1:2:2*nn1)  = N1;
        Nm1(2,2:2:2*nn1)  = N1;
        
     
        B2(1,1:2:2*nn2)  = (-beta/3*y^3) *dN2dx2;
        B2(1,2:2:2*nn2)  = (-y+beta/3*y^3)*dN2dx;
        B2(2,1:2:2*nn2)  = (1-beta*y^2)  *dN2dx;
        B2(2,2:2:2*nn2)  = (1-beta*y^2)  *N2;
                       
        Nm2(1,1:2:2*nn2)  =  (-beta/3*y^3) *dN2dx;
        Nm2(1,2:2:2*nn2)  =  (-y+beta/3*y^3)*N2;
        Nm2(2,1:2:2*nn2)  = N2;
        
        
        yPt=N1*pts1;
        
        % COMPUTE ELEMENT STRAIN AND STRESS AT STRESS POINT
        strain=B1*U(sctrB1);
        stress=C*strain;
        sigmaSolid(e,1)    = stress(1);
        sigmaSolid(e,2)    = stress(3);
        sigmaBeam (e,:)    = Cbeam*B2*U(sctrB2);
        sigmaRef(e,1) = 1000/I0*(L-yPt(1))*yPt(2);
        sigmaRef(e,2) = -1000/2/I0*(c^2-yPt(2)^2);
        xcoord(e,:)        = yPt;
    
        dispRef(e,1) = 1000*yPt(2)/(6*E0*I0)*((6*L-3*yPt(1))*yPt(1)+(2+nu0)*(yPt(2)^2-D^2/4));
        dispRef(e,2) = -1000/(6*E0*I0)*(3*nu0*yPt(2)^2*(L-yPt(1))+(4+5*nu0)*D^2*yPt(1)/4+...
            (3*L-yPt(1))*yPt(1)^2);
        
        dispSolid(e,:)  = Nm1 * U(sctrB1);
        dispBeam(e,:)   = Nm2 * U(sctrB2);

        e = e + 1;
    end  % of quadrature loop
end

colordef white
figure,set (gcf,'Color','w')
set(gca,'FontSize',14)
hold on
plot(xcoord(:,2),sigmaSolid(:,2),'o','MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',6.5);
%plot(xcoord(:,2),sigmaBeam(:,2),'s','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',6.5);
plot(xcoord(:,2),sigmaRef(:,2),'*','MarkerEdgeColor','k','MarkerFaceColor','cy','MarkerSize',6.5);
% plot(xcoord(:,2),sigmaRef(:,2),'k-','LineWidth',1.4);
% plot(xcoord(:,2),sigma(:,2),'s','MarkerEdgeColor','k',...
%     'MarkerFaceColor','g',...
%     'MarkerSize',6.5);
h=legend('solid','beam');
xlabel('y')
ylabel('sigma_{xy} at x=L/2')
grid on

colordef white
figure,set (gcf,'Color','w')
set(gca,'FontSize',14)
hold on
plot(xcoord(:,2),dispSolid(:,2),'o','MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',6.5);
%plot(xcoord(:,2),dispBeam (:,2),'s','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',6.5);
%plot(xcoord(:,2),dispRef  (:,2),'*','MarkerEdgeColor','k','MarkerFaceColor','cy','MarkerSize',6.5);
h=legend('solid','beam','reference');
xlabel('y')
ylabel('sigma_{xy} at x=L/2')
grid on
