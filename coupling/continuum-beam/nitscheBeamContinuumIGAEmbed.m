% This file implements the Nitsche method to join two mechanical models:
% a 2D continuum model and a Timoshenko beam model.
% Embeded mesh method in which the 2D continuum model is embedded on a
% background beam mesh.
%
% Discretisation: B-spline elements.
%
% Problem: Timoshenko beam in bending.
%
% Vinh Phu Nguyen
% Cardiff University, UK
% 12 July 2013
% 8 August 2013: remedy for singular matrix when interface is close to
% the beam cut element boundary.


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


%% MATERIAL PROPERTIES-----------------------------------------------------
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
shearFactor = 5/6;

plotMesh  = 1;

stressState = 'PLANE_STRESS';
% COMPUTE ELASTICITY MATRIX
C = elasticityMatrix(E0,nu0,stressState);

% penalty parameter

alpha = 1e10;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GENERATE FINITE ELEMENT MESH
%

disp([num2str(toc),'   GENERATING MESH'])

% MESH PROPERTIES

L1 = 29.97;
% meshing for domain1------------------------------------------------------

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

solid1 = nrbdegelev(solid1,[0 0]);

uKnot     = cell2mat(solid1.knots(1));
vKnot     = cell2mat(solid1.knots(2));

refineCountX = 2;
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

refineCountX = 5;
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

% domain 2-----------------------------------------------------------------

controlPts          = zeros(4,2);

controlPts(1:2,1) = [0;0];
controlPts(1:2,2) = [L;0];

controlPts(4,:)   = 1;

uKnot = [0 0 1 1];

solid2 = nrbmak(controlPts,uKnot);

% evaluate order

solid2 = nrbdegelev(solid2,3);

uKnot     = solid2.knots;

refineCountX = 3;
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
    %plot(mesh1.controlPts(:,1),mesh1.controlPts(:,2),'o','MarkerEdgeColor','k',...
    %    'MarkerFaceColor','g','MarkerSize',6.5);
    nrbkntplot(solid2);
    plot(mesh2.controlPts(:,1),mesh2.controlPts(:,2),'o','MarkerEdgeColor','k',...
        'MarkerFaceColor','w','MarkerSize',6.5);
    hold on
    axis off
    axis([0 L -c c])
end

% boundary mesh for domain 1
bndMesh1 = [];
for i=1:mesh1.noElemsV
    bndMesh1 = [bndMesh1; mesh1.noElemsU*i ];
end

% boundary edges

bndNodes  = find(abs(mesh1.controlPts(:,1)-L1)<1e-12);
bndEdge1  = zeros(mesh1.noElemsV,mesh1.q+1);

for i=1:mesh1.noElemsV
    bndEdge1(i,:) = bndNodes(i:i+mesh1.q);
end

%%

xi = (L1-0)/L; % parameter space of L1

knotSpanIndex = FindSpan(solid2.number-1,solid2.order-1,xi,solid2.knots);
cutElem = knotSpanIndex-solid2.order+2;

elem2State = zeros(size(mesh2.element,1),1);
elem2State(:) = 0;
id = find(1:mesh2.elemCount >= cutElem);
elem2State(id) = 2;
elem2State(cutElem) = 1;

% find dofs of beam that are complete covered by the solid mesh
% these dofs are named inactive dofs and will be constrained.

inactiveNodes1 = [];

for i=1:size(mesh2.element,1)
    sctr = mesh2.element(i,:);
    if elem2State(i) == 0
        inactiveNodes1 = [inactiveNodes1;sctr];
    end
end

inactiveNodes1 = unique(inactiveNodes1);
inactiveNodes2 = mesh2.element(cutElem,:);
inactiveNodes2 = inactiveNodes2(:);

%%%% !!!!
% if the cut element is almost covered by the solid then
% mark the first control point of cut element as inactive 
% because support domain of this control point is covered by the solid.
inactiveNodes2(1)=[];

inactiveBeamNodes = setdiff(inactiveNodes1,inactiveNodes2);
inactiveBeamNodes = inactiveBeamNodes + size(mesh1.controlPts,1);

%% Gauss points on coupling interface--------------------------------------
ngp = 9;
[W1,Q1]=quadrature( ngp, 'GAUSS', 1 ); % two point quadrature

GP1 = [];
GP2 = [];

for e=1:mesh1.noElemsV
    sctrEdge = bndEdge1(e,:);
    sctr1    = mesh1.globElems(bndMesh1(e),:);
    pts1     = mesh1.controlPts(sctr1,:);
    xiE1     = mesh1.rangeU(end,:); % [xi_i,xi_i+1]
    etE1     = mesh1.rangeV(e,:);   % [et_j,eta_j+1]
    xiE2     = mesh2.range(1,:);   % [et_j,eta_j+1]
    for q=1:size(W1)
        pt=Q1(q,:);
        wt=W1(q);
        Xi      = 0.5 * ( ( etE1(2) - etE1(1) ) * pt + etE1(2) + etE1(1));
        J2      = 0.5 * ( etE1(2) - etE1(1) );
        [N dNdxi] = NURBS1DBasisDers(Xi,mesh1.q,mesh1.vKnot,mesh1.weights);
        J0=dNdxi*mesh1.controlPts(sctrEdge,:); % also the tangent to Gamma_*^e
        detJ0=norm(J0);
        J0 = J0/detJ0;
        
        x  = N*mesh1.controlPts(sctrEdge,:);
        X1 = global2LocalMapNURBS2D(x,pts1,xiE1,etE1,mesh1.p,mesh1.q, ...
            mesh1.uKnot, mesh1.vKnot, mesh1.weights);
        X2  = x(1)/L;
        GP1 = [GP1;X1 wt*detJ0*J2 -J0(2) J0(1)];
        GP2 = [GP2;X2 x(2)];
    end
end

% DEFINE BOUNDARIES

% GET NODES ON DISPLACEMENT BOUNDARY
%      Here we get the nodes on the essential boundaries

fixedNode = find(mesh1.controlPts(:,1)==0);
rightNode = find(mesh2.controlPts(:,1)==L);
midNode1  = find(mesh1.controlPts(:,2)==0);

%% least square for boujndary conditions

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

noxC   = 8;

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

fixedNode = [fixedNode;inactiveBeamNodes];

uFixed = [uFixed;zeros(length(inactiveBeamNodes),1)];
vFixed = [vFixed;zeros(length(inactiveBeamNodes),1)];

numdofs = size(mesh1.controlPts,1) + size(mesh2.controlPts,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE SYSTEM DATA STRUCTURES

disp([num2str(toc),' INITIALIZING DATA STRUCTURES'])

f=zeros(2*numdofs,1);          % external load vector
K=zeros(2*numdofs,2*numdofs); % stiffness matrix

% ******************************************************************************
% ***                          P R O C E S S I N G                           ***
% ******************************************************************************
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%% COMPUTE STIFFNESS MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp([num2str(toc),'   COMPUTING STIFFNESS MATRIX'])

%% for domain 1------------------------------------------------------------

ngpv = mesh1.p+1;

[W,Q]=quadrature( ngpv, 'GAUSS', 2 ); % 2x2 Gaussian quadrature

for e=1:mesh1.elemCount
    idu    = mesh1.index(e,1);
    idv    = mesh1.index(e,2);
    xiE    = mesh1.rangeU(idu,:); % [xi_i,xi_i+1]
    etaE   = mesh1.rangeV(idv,:); % [eta_j,eta_j+1]
    
    sctr   = mesh1.globElems(e,:);         %  element scatter vector
    nn     = length(sctr);
    
    sctrB(1,1:2:2*nn) = 2*sctr-1;
    sctrB(1,2:2:2*nn) = 2*sctr  ;
    
    pts    = mesh1.controlPts(sctr,:);
    
    % loop over Gauss points
    
    for gp=1:size(W,1)
        pt      = Q(gp,:);
        wt      = W(gp);
        Xi      = parent2ParametricSpace(xiE,pt(1));
        Eta     = parent2ParametricSpace(etaE,pt(2));
        J2      = jacobianPaPaMapping(xiE,etaE);
        
        % compute derivatives of shape functions
        [R dRdxi dRdeta] = NURBS2DBasisDers([Xi; Eta],mesh1.p,mesh1.q,mesh1.uKnot,mesh1.vKnot,mesh1.weights');
        
        % compute the jacobian of physical and parameter domain mapping
        % then the derivative w.r.t spatial physical coordinates
        
        jacob      =  pts' * [dRdxi' dRdeta'];
        J1         = det(jacob);
        dRdx       = [dRdxi' dRdeta'] * inv(jacob);
        J          = J1 * J2;
        B           = getBmatrix2D(dRdx');
        K(sctrB,sctrB) = K(sctrB,sctrB) + B' * C * B * J * wt;
    end
end

%% for domain 2 (beam)-----------------------------------------------------
clear sctrB;
noGPs=solid2.order;
[W,Q]=quadrature(  noGPs, 'GAUSS', 1 ); % 2 point quadrature
xi = (L1-0)/L; % parameter space of L1
for e=1:mesh2.elemCount
   if elem2State(e) == 0
        continue
   end
   
   if e==cutElem
       [W,Q]=quadrature(  10, 'GAUSS', 1 ); % 2 point quadrature
   end
   
   xiE   = mesh2.range(e,:); % [xi_i,xi_i+1]
   conn  = mesh2.element(e,:);
   noFns = length(conn);
   conng = conn + size(mesh1.controlPts,1);
   sctrB(1:2:2*noFns) = 2*conng-1;
   sctrB(2:2:2*noFns) = 2*conng-0;
   
   pts   = mesh2.controlPts(conn,1:2);
   
   % loop over Gauss points 
    for gp=1:size(W,1)                        
      pt      = Q(gp,:);                          
      wt      = W(gp);                 
      % coord in parameter space  
      Xi      = 0.5 * (( xiE(2) - xiE(1) ) * pt + xiE(2) + xiE(1));
      J2      = 0.5 * ( xiE(2) - xiE(1) ); 
      
      if ( e == cutElem && Xi < xi )
          continue
      end
      
      % shape functions and first derivatives 
            
      [R dRdxi] = NURBS1DBasisDers(Xi,solid2.order-1,uKnot,solid2.coefs(4,:));
      
      % compute the jacobian of physical and parameter domain mapping
      % then the derivative w.r.t spatial physical coordinates
      
      dxdxi   = dRdxi *pts;      
     
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

% interface integrals

Cbeam  = [E0 0;0 shearFactor*G];
Csolid = C([1 3],:);

for i=1:mesh1.noElemsV                     % start of element loop
    e1     = bndMesh1(i);
    sctr1  = mesh1.globElems(e1,:);
    sctr2  = mesh2.element(cutElem,:);
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
        
        B2(1,2:2:2*nn2)  = -y*dN2dx;
        B2(2,1:2:2*nn2)  = dN2dx;
        B2(2,2:2:2*nn2)  = -N2;
        
        Nm1(1,1:2:2*nn1)  = N1';
        Nm1(2,2:2:2*nn1)  = N1';
        
        Nm2(1,2:2:2*nn2)  = -y*N2';
        Nm2(2,1:2:2*nn2)  = N2';
        
        % COMPUTE ELEMENT STIFFNESS AT QUADRATURE POINT
        
        Kp11 = Kp11 + alpha*(Nm1'*Nm1)*wt1;
        Kp12 = Kp12 + alpha*(Nm1'*Nm2)*wt1;
        Kp22 = Kp22 + alpha*(Nm2'*Nm2)*wt1;
        
        Kd11 = Kd11 + 0.5 * Nm1'* n * Csolid * B1 *wt1;
        Kd12 = Kd12 + 0.5 * Nm1'* n * Cbeam  * B2 *wt1;
        Kd21 = Kd21 + 0.5 * Nm2'* n * Csolid * B1 *wt1;
        Kd22 = Kd22 + 0.5 * Nm2'* n * Cbeam  * B2 *wt1;
        
    end  % of quadrature loop
    
    K(sctrB1,sctrB1)  = K(sctrB1,sctrB1)  - Kd11 - Kd11' + Kp11;
    K(sctrB1,sctrB2)  = K(sctrB1,sctrB2)  - Kd12 + Kd21' - Kp12;
    K(sctrB2,sctrB1)  = K(sctrB2,sctrB1)  + Kd21 - Kd12' - Kp12';
    K(sctrB2,sctrB2)  = K(sctrB2,sctrB2)  + Kd22 + Kd22' + Kp22;
end


f(2*numdofs-1) = P;

%%%%%%%%%%%%%%%%%%% END OF STIFFNESS MATRIX COMPUTATION %%%%%%%%%%%%%%%%%%%


%% APPLY ESSENTIAL BOUNDARY CONDITIONS
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
f(vdofs)=bcwt*speye(length(udofs))*vFixed;

%% SOLVE SYSTEM
disp([num2str(toc),'   SOLVING SYSTEM'])
U=K\f;

%******************************************************************************
%*** POST - PROCESSING ***
%***************************************************

numnode1 = size(mesh1.controlPts,1);
numnode2 = size(controlPts,1);

Ux = U(1:2:2*numdofs);
Uy = U(2:2:2*numdofs);

Ux1 = Ux(1:numnode1);
Ux2 = Ux(numnode1+1:end);

Uy1 = Uy(1:numnode1);
Uy2 = Uy(numnode1+1:end);

%%

Uxm1  = U(2*midNode1-1);
Uym1  = U(2*midNode1);

Ux2(1:length(inactiveBeamNodes))=[];

xx = mesh1.controlPts(midNode1,1);
u  = Uym1;

xx2 = [mesh2.controlPts(length(inactiveBeamNodes)+1:end,1)];
u2  = [Ux2];

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
plot(xx,u,'o','MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',6.5);
plot(xx2,Ux2,'s','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',6.5);
h=legend('exact','continuum','beam');
xlabel('x')
ylabel('w')
grid on

%%
% factor=100;
% ctrpts1 = reshape(solid1.coefs,4,[]).';
% ctrpts1(:,[1 2]) = ctrpts1(:,[1 2]) + factor*[Ux1 Uy1];
% ctrpts1 = reshape(ctrpts1.', 4, solid1.number(1), solid1.number(2));
% solid1 = nrbmak(ctrpts1,{solid1.knots{1} solid1.knots{2}});
% 
% ctrpts2 = reshape(solid2.coefs,4,[])';
% ctrpts2(:,[1 2]) = ctrpts2(:,[1 2]) + factor*[Uy2 Ux2];
% ctrpts2 = reshape(ctrpts2.', 4, solid2.number);
% solid2  = nrbmak(ctrpts2,solid2.knots);
% 
% figure
% hold on
% nrbkntplot(solid1);
% nrbkntplot(solid2);
% view(0,90);

%%
clear sctrB;
aa=floor(mesh1.noElemsV/2);
elems1 = mesh1.noElemsU*aa+1:mesh1.noElemsU*(aa+1);
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
        pt      = [0 0];
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
plot(xcoord(:,1),sigmaRef(:,1),'k-','LineWidth',1.4);
plot(xcoord(:,1),sigma(:,1),'o','MarkerEdgeColor','k',...
    'MarkerFaceColor','g',...
    'MarkerSize',6.5);
plot(xcoord(:,1),sigmaRef(:,2),'k-','LineWidth',1.4);
plot(xcoord(:,1),sigma(:,2),'s','MarkerEdgeColor','k',...
    'MarkerFaceColor','g',...
    'MarkerSize',6.5);
h=legend('sigmaxx-exact','sigmaxx-coupling','sigmaxy-exact','sigmaxy-coupling');
xlabel('x')
ylabel('sigma_{xy} at y=L/2')
grid on
    
%%

aa=mesh1.noElemsU;
elems1 = aa:mesh1.noElemsU:aa+mesh1.noElemsU*(mesh1.noElemsV-1);

sigma      = zeros(length(elems1),2);
sigmaRef   = zeros(length(elems1),2);
xcoord     = zeros(length(elems1),2);

e = 1;
for i=1:mesh1.noElemsV                     % start of element loop
    e1     = bndMesh1(i);
    sctr1  = mesh1.globElems(e1,:);        
    nn1    = length(sctr1);    
    sctrB1(1:2:2*nn1) = 2*sctr1-1;
    sctrB1(2:2:2*nn1) = 2*sctr1-0;
    
    pts1 = mesh1.controlPts(sctr1,:);

    for q=ngp*(i-1)+1:ngp*(i-1)+ngp        
        pt1 = GP1(q,1:2);
        wt1 = GP1(q,3);
        [N1, dN1dxi, dN1deta] = NURBS2DBasisDers(pt1,mesh1.p,mesh1.q,mesh1.uKnot,mesh1.vKnot,mesh1.weights');
                
        J1        = pts1' * [dN1dxi' dN1deta'];                        
        dN1dx     = [dN1dxi' dN1deta'] * inv(J1);                
        
        B1(1,1:2:2*nn1)  = dN1dx(:,1)';
        B1(2,2:2:2*nn1)  = dN1dx(:,2)';
        B1(3,1:2:2*nn1)  = dN1dx(:,2)';
        B1(3,2:2:2*nn1)  = dN1dx(:,1)';   
        
        yPt=N1*pts1;
        
        % COMPUTE ELEMENT STRAIN AND STRESS AT STRESS POINT
        strain=B1*U(sctrB1);
        stress=C*strain;
        sigma(e,1)    = stress(1);
        sigma(e,2)    = stress(3);
        sigmaRef(e,1) = 1000/I0*(L-yPt(1))*yPt(2);
        sigmaRef(e,2) = -1000/2/I0*(c^2-yPt(2)^2);
        xcoord(e,:)     = yPt;
    
        e = e + 1;
    end  % of quadrature loop
end

colordef white
figure,set (gcf,'Color','w')
set(gca,'FontSize',14)
hold on
plot(xcoord(:,2),sigmaRef(:,2),'k-','LineWidth',1.4);
plot(xcoord(:,2),sigma(:,2),'o','MarkerEdgeColor','k',...
    'MarkerFaceColor','g',...
    'MarkerSize',6.5);
% plot(xcoord(:,2),sigmaRef(:,2),'k-','LineWidth',1.4);
% plot(xcoord(:,2),sigma(:,2),'s','MarkerEdgeColor','k',...
%     'MarkerFaceColor','g',...
%     'MarkerSize',6.5);
h=legend('sigmaxy-exact','sigmaxy-coupling');
xlabel('y')
ylabel('sigma_{xy} at x=L/2')
grid on

%%


% Here we plot the stresses and displacements of the solution. As with the
% mesh generation section we don?t go into too much detail - use help
% ?function name? to get more details.
disp([num2str(toc),'   POST-PROCESSING'])

material.stiffMat=C;
displacement = [U(1:2:end) U(2:2:end)];

matMap1 = ones(mesh1.elemCount,1);

vMesh1=buildVisualizationMesh2D(solid1);

meshes{1}  = mesh1;
vmeshes{1} = vMesh1;

meshData.mesh    = meshes;
meshData.vmesh   = vmeshes;
% meshData.noElems = meshes{1}.noElems;
% meshData.noPts   = meshes{1}.noPts;
meshData.matMap{1}=matMap1;
materials{1}      = material;

vtuFileName = 'iga-nitsche-2D-beam';
for ip=1:1
    vtuFile = strcat(vtuFileName,'-mesh',num2str(ip));
    %figure; hold on;
    ok      = plotStress2DForPatch(meshData,ip,vtuFile,displacement,materials);
end


