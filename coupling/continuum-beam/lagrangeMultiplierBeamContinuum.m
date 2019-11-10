% This file implements the Lagrange multiplier method to join two mechanical models:
% a 2D continuum model and a Timoshenko beam model.
% This serves exact solution to which Nitsche's solution should be compared
% with.
%
% Discretisation: standard Lagrange finite elements.
%
% Problem: Timoshenko beam in bending.
%
% Vinh Phu Nguyen
% Cardiff University, UK
% 4 November 2013

addpath ../../fem_util/
addpath ../../gmshFiles/
addpath ../../post-processing/
addpath ../../fem-functions/
addpath ../../analytical-solutions/


clear all
clc
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
q  = 0; % uniformly distributed loads
P  = -1000*t^3/12/I0; % end force
G  = E0/2/(1+nu0);
k  = E0*I0;
shearFactor = 5/6;

plotMesh  = 1;

stressState = 'PLANE_STRESS';
% COMPUTE ELASTICITY MATRIX
C = elasticityMatrix(E0,nu0,stressState);

noGPs = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE FINITE ELEMENT MESH
%

disp([num2str(toc),'   GENERATING MESH'])

% MESH PROPERTIES

elemType  = 'Q4'; % the element type used in the solid domain;
elemTypeB = 'L2'; % the element type used in the solid domain (boundary shape);

%% domain 1 ---------------------------------------------------------------
numx1     = 40;
numy1     = 9; % should be an odd number so that stresses at the midline 
               % can be easily obtained (see post-processing for details)

L1 = 24;
% meshing for domain1
nnx=numx1+1;
nny=numy1+1;
node1=square_node_array([0 -c],[L1 -c],[L1 c],[0 c],nnx,nny);
inc_u=1;
inc_v=nnx;
node_pattern=[ 1 2 nnx+2 nnx+1 ];
element1 =make_elem(node_pattern,numx1,numy1,inc_u,inc_v);

%% domain 2 ---------------------------------------------------------------

numnode = 16;
node2 = zeros(numnode,2);
node2(:,1) = linspace(L1,L,numnode);

element2  = zeros(numnode-1,2);

for i=1:numnode-1
    element2(i,:) = i:i+1;
end


% PLOT MESH
if ( plotMesh )
    clf
    plot_mesh(node1,element1,elemType,'g.-',2.4);
    plot_mesh(node2,element2,elemTypeB,'r.-',1.7);
    hold on
    axis off
    axis([0 L -c c])
end

% boundary mesh for domain 1
bndMesh1 = [];
for i=1:numy1
    bndMesh1 = [bndMesh1; numx1*i ];
end

% boundary edges

bndNodes  = find(abs(node1(:,1)-L1)<1e-14);
bndEdge1  = zeros(numy1,2);

for i=1:numy1
    bndEdge1(i,:) = bndNodes(i:i+1);
end

ngp=2;
[W1,Q1]=quadrature( ngp, 'GAUSS', 1 ); % two point quadrature

GP1 = [];
GP2 = [];

for e=1:numy1
    sctrEdge=bndEdge1(e,:);
    sctr1    = element1(bndMesh1(e),:);
    pts1     = node1(sctr1,:);
    for q=1:size(W1)
        pt=Q1(q,:);
        wt=W1(q);
        [N,dNdxi]=lagrange_basis('L2',pt);  % element shape functions
        J0=dNdxi'*node1(sctrEdge,:);
        detJ0=norm(J0);
        J0 = J0/detJ0;
        
        x=N'*node1(sctrEdge,:);
        X1 = global2LocalMap(x,pts1,elemType);
        GP1 = [GP1;X1 wt*detJ0 -J0(2) J0(1) pt];
        GP2 = [GP2;-1 x(2)];
    end
end

% DEFINE BOUNDARIES

edgeElemType='L2';

% GET NODES ON DISPLACEMENT BOUNDARY
%      Here we get the nodes on the essential boundaries

fixedNode = find(node1(:,1)==0);
rightNode = find(node2(:,1)==L);
midNode1  = find(node1(:,2)==0);

uFixed=zeros(1,length(fixedNode))';  % a vector of u_x for the nodes
vFixed=zeros(1,length(fixedNode))';

for i=1:length(fixedNode)
    inode = fixedNode(i);
    pts   = node1(inode,:);
    ux    = 1000*pts(2)/(6*E0*I0)*(2+nu0)*(pts(2)^2-c^2);
    uy    = -1000/(2*E0*I0)*nu0*pts(2)^2*L;
    uFixed(i)=ux;
    vFixed(i)=uy;
end

udofs=2*fixedNode-1;           % global indecies of the fixed x displacements
vdofs=2*fixedNode;             % global indecies of the fixed y displacements

numdofs = size(node1,1) + size(node2,1);

inodes  = 1:numdofs;

fiNodes=[fixedNode;numdofs];
frNodes=setdiff(1:numdofs,fiNodes);

activeDofs=[];
for i=1:length(frNodes)
    in = frNodes(i);
    activeDofs=[activeDofs; 2*in-1;2*in];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE SYSTEM DATA STRUCTURES

disp([num2str(toc),' INITIALIZING DATA STRUCTURES'])

f = zeros(2*numdofs,1);          % external load vector
K = zeros(2*numdofs,2*numdofs);  % stiffness matrix


% ******************************************************************************
% ***                          P R O C E S S I N G                           ***
% ******************************************************************************
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%% COMPUTE STIFFNESS MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp([num2str(toc),'   COMPUTING STIFFNESS MATRIX'])

[W,Q]=quadrature( 2, 'GAUSS', 2 ); % 2x2 Gaussian quadrature

%% for domain 1------------------------------------------------------------

for e=1:size(element1,1)                          % start of element loop
    sctr = element1(e,:);           % element scatter vector
    nn   = length(sctr);
    sctrB(1:2:2*nn) = 2*sctr-1;
    sctrB(2:2:2*nn) = 2*sctr-0;
    
    for q=1:size(W,1)
        pt=Q(q,:);
        wt=W(q);
        [N,dNdxi]=lagrange_basis(elemType,pt); % element shape functions
        J0=node1(sctr,:)'*dNdxi;                % element Jacobian matrix
        invJ0=inv(J0);
        dNdx=dNdxi*invJ0;
        
        B(1,1:2:2*nn)  = dNdx(:,1)';
        B(2,2:2:2*nn)  = dNdx(:,2)';
        B(3,1:2:2*nn)  = dNdx(:,2)';
        B(3,2:2:2*nn)  = dNdx(:,1)';
        
        % COMPUTE ELEMENT STIFFNESS AT QUADRATURE POINT
        K(sctrB,sctrB)=K(sctrB,sctrB)+B'*C*B*W(q)*det(J0);
    end  % of quadrature loop
end

Ks = K(1:2*size(node1,1),1:2*size(node1,1));

%% for domain 2 (beam)-----------------------------------------------------
clear sctrB;

[W,Q]=quadrature(  noGPs, 'GAUSS', 1 ); % 2 point quadrature

for e=1:size(element2,1)
    conn  = element2(e,:);
    noFns = length(conn);
    conng = conn + size(node1,1);
    sctrB(1:2:2*noFns) = 2*conng-1;
    sctrB(2:2:2*noFns) = 2*conng-0;
    pts   = node2(conn,:);
    
    % loop over Gauss points
    for gp=1:size(W,1)
        pt      = Q(gp,:);
        wt      = W(gp);
        [N,dNdxi]=lagrange_basis(elemTypeB,pt);  % element shape functions
        J0=dNdxi'*pts;
        detJ0=norm(J0);
        dNdx   = (1/detJ0)*dNdxi;
        
        Bb(2:2:noFns*2)   = dNdx;
        Bs(1:2:noFns*2)   = dNdx;
        Bs(2:2:noFns*2)   = -N;
        
        % compute elementary stiffness matrix
        K(sctrB,sctrB) = K(sctrB,sctrB) + ...
            (k * Bb' * Bb + shearFactor*G*b*t*Bs'*Bs)* detJ0 * wt;
    end
end

%% Lagrange multiplier method
% there are three unknowns: solid displacements (S), beam displacements (B)
% and Lagrange multipliers (L).
% In the following, it is assumed that L is approximated as S.

num_disp_nodes = 2*length(bndNodes);
G = zeros(2*numdofs,num_disp_nodes); % Lagrange matrix

beta=0;
for i=1:numy1                     % start of element loop
    e1     = bndMesh1(i);
    sctrS  = element1(e1,:);            % connectivity of solid element
    sctrB  = element2(1,:);             % connectivity of beam element
    sctrL  = i:i+1;                     % connectivity of Lagrange multiplier element
    sctrBn = sctrB + size(node1,1);
    nn1    = length(sctrS);
    nn2    = length(sctrB);
    nn3    = length(sctrL);
    
    sctrBS(1:2:2*nn1) = 2*sctrS-1;
    sctrBS(2:2:2*nn1) = 2*sctrS-0;
    sctrBB(1:2:2*nn2) = 2*sctrBn-1;
    sctrBB(2:2:2*nn2) = 2*sctrBn-0;
    sctrBL(1:2:2*nn3) = 2*sctrL-1;
    sctrBL(2:2:2*nn3) = 2*sctrL-0;
    
    pts1 = node1(sctrS,:);
    pts2 = node2(sctrB,:);
        
    for q=ngp*(i-1)+1:ngp*(i-1)+ngp  % loop over Gauss points 
        
        pt1=GP1(q,1:2);
        wt1=GP1(q,3);      
        pt2=GP2(q,1);
        y  =GP2(q,2);
        pt3=GP1(q,6);
        
        [N1,dN1dxi] = lagrange_basis(elemType, pt1);
        [N2,dN2dxi] = lagrange_basis(elemTypeB,pt2);
        [NL,dNLdxi] = lagrange_basis('L2'     ,pt3);
        
        J1 = pts1'*dN1dxi;
        J2 = pts2'*dN2dxi;
        detJ0=norm(J2);
        
        dN1dx   = dN1dxi*inv(J1);
        dN2dx   = (1/detJ0)*dN2dxi;
                
        NmS(1,1:2:2*nn1)  = N1';
        NmS(2,2:2:2*nn1)  = N1';
                               
        NmB(1,1:2:2*nn2)  =  (-beta/3*y^3) *dN2dx;
        NmB(1,2:2:2*nn2)  =  (-y+beta/3*y^3)*N2;
        NmB(2,1:2:2*nn2)  = N2;
        
        NmL(1,1:2:2*nn3)  = NL';
        NmL(2,2:2:2*nn3)  = NL';
                                               
        G(sctrBS,sctrBL) = G(sctrBS,sctrBL) + NmS'*NmL * wt1;
        G(sctrBB,sctrBL) = G(sctrBB,sctrBL) - NmB'*NmL * wt1;        
    end  % of quadrature loop
end


f(2*numdofs-1) = P;

qk   = zeros(1,num_disp_nodes);
fnew = [f;qk'];                                % f = {f;qk}
Knew = ([K G; G' zeros(num_disp_nodes)]);    % m = [K GG;GG' 0]

%%%%%%%%%%%%%%%%%%% END OF STIFFNESS MATRIX COMPUTATION %%%%%%%%%%%%%%%%%%%


% APPLY ESSENTIAL BOUNDARY CONDITIONS
disp([num2str(toc),'   APPLYING BOUNDARY CONDITIONS'])
bcwt=mean(diag(Knew)); % a measure of the average size of an element in K
% used to keep the conditioning of the K matrix

fnew=fnew-Knew(:,udofs)*uFixed;  % modify the force vector
fnew=fnew-Knew(:,vdofs)*vFixed;
Knew(udofs,:)=0;
Knew(vdofs,:)=0;
Knew(:,udofs)=0;
Knew(:,vdofs)=0;
Knew(udofs,udofs)=bcwt*speye(length(udofs)); % put ones*bcwt on the diagonal
Knew(vdofs,vdofs)=bcwt*speye(length(vdofs));
fnew(udofs)=bcwt*speye(length(udofs))*uFixed;
fnew(vdofs)=bcwt*speye(length(udofs))*vFixed;

% eigenvalues of K to check positive definiteness

% eigens = eig(K);
% max(eigens)

% SOLVE SYSTEM
disp([num2str(toc),'   SOLVING SYSTEM'])
solution=Knew\fnew;

U=solution(1:2*numdofs);

%******************************************************************************
%*** POST - PROCESSING ***
%***************************************************

numnode1 = size(node1,1);
numnode2 = size(node2,1);

Ux = U(1:2:2*numdofs);
Uy = U(2:2:2*numdofs);

Ux1 = Ux(1:numnode1);
Ux2 = Ux(numnode1+1:end);

Uy1 = Uy(1:numnode1);
Uy2 = Uy(numnode1+1:end);


disp([num2str(toc),'   POST-PROCESSING'])


fn=1;

%% exact solution at mesh1 nodes

exactDisp = zeros(size(node1,1),2);
for i=1:size(node1,1)
    x = node1(i,1); 
    y = node1(i,2);
    
    % exact displacements    
    ux = 1000*y/(6*E0*I0)*((6*L-3*x)*x+(2+nu0)*(y^2-t^2/4));
    uy = -1000/(6*E0*I0)*(3*nu0*y^2*(L-x)+(4+5*nu0)*t^2*x/4+(3*L-x)*x^2);
    exactDisp(i,:) = [ux uy];
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT DEFORMED DISPLACEMENT PLOT
%Ux = U(xs);
%Uy = U(ys);

% scaleFact=300.;
% colordef white
% figure
% clf
% hold on
% %plot_field(node1+scaleFact*[Ux1 Uy1],element1,elemType,Ux1);
% %plot_field(node2+scaleFact*[Ux2 Uy2],element2,elemType,Ux2);
% plot_mesh(node1+scaleFact*[Ux1 Uy1],element1,elemType,'w.-',1);
% plot_mesh(node2+scaleFact*[zeros(numnode2,1) Ux2],element2,elemTypeB,'w.-',1);
% plot_mesh(node1+scaleFact*exactDisp,element1,elemType,'k.-',1);
% axis off
% %title('DEFORMED DISPLACEMENT IN Y-DIRECTION')
% 

%% Comapre numerical displacement to exact value
[W,Q]=quadrature(  10, 'GAUSS', 1 ); % 2 point quadrature

Uxm1  = U(2*midNode1-1);
Uym1  = U(2*midNode1);
xx    = node1(midNode1,1);

% if numy1=odd number, thre is no nodes along the midline
% we have to compute the displacement from nodal values
if isempty(midNode1)
    no    = floor(numy1/2);
    elems = numx1*(no)+1:numx1*(no+1);
    
    for e=1:length(elems)
        ie = elems(e);
        sctr=element1(ie,:);
        nn   = length(sctr);
        sctrx = 2*sctr-1;
        sctry = 2*sctr-0;
        nn=length(sctr);
        
        pt=[0 0];
        [N,dNdxi]=lagrange_basis(elemType,pt);

        yPt=N'*node1(sctr,:);
        Uxm1 = [Uxm1; N'*U(sctrx)];
        Uym1 = [Uym1; N'*U(sctry)];
        xx   = [xx;yPt(1)];
    end   % of element loop
end


u  = Uym1;

xx = [xx; node2(:,1)];
u  = [u;  Ux2];

% for e=1:size(element2,1)
%    conn  = element2(e,:);
%    noFns = length(conn);
%    conng = conn + size(node1,1);
%    sctrx = 2*conng-1;
%    sctry = 2*conng-0;
%    pts   = node2(conn,:);
% 
%    % loop over Gauss points
%     for gp=1:size(W,1)
%       pt      = Q(gp,:);
%       wt      = W(gp);
%       [N,dNdxi]=lagrange_basis(elemTypeB,pt);  % element shape functions
%       J0=dNdxi'*pts;
%       xx = [xx; N'*pts(:,1)];
%       u  = [u;  N'*U(sctrx)];
%     end
% end


% xx=[node2(:,1)];
% u=[Ux2];


% exact solution u =
% y= 0;
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
% plot(uCont(:,1),uCont(:,2),'b--','LineWidth',1.4);
% plot(uCoup(:,1),uCoup(:,2),'cy.-','LineWidth',1.4);
h=legend('exact','coupling');
xlabel('x')
ylabel('w')
grid on
%axis([0 5 -0.55 0])


% aa=40;
% elems = aa:numx1:aa+numx1*(numy1-1);
% 
% sigma      = zeros(length(elems),2);
% disp       = zeros(length(elems),2);
% sigmaRef   = zeros(length(elems),2);
% xcoord     = zeros(length(elems),2);
% 
% for e=1:length(elems)
%     ie = elems(e);
%     sctr=element1(ie,:);
%     nn   = length(sctr);
%     sctrB(1:2:2*nn) = 2*sctr-1;
%     sctrB(2:2:2*nn) = 2*sctr-0;
%     nn=length(sctr);
%     
%     pt=[0 0];
%     [N,dNdxi]=lagrange_basis(elemType,pt);
%     J0=node1(sctr,:)'*dNdxi;
%     invJ0=inv(J0);
%     dNdx=dNdxi*invJ0;
%     yPt=N'*node1(sctr,:);
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % COMPUTE B MATRIX
%     B(1,1:2:2*nn)  = dNdx(:,1)';
%     B(2,2:2:2*nn)  = dNdx(:,2)';
%     B(3,1:2:2*nn)  = dNdx(:,2)';
%     B(3,2:2:2*nn)  = dNdx(:,1)';
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % COMPUTE ELEMENT STRAIN AND STRESS AT STRESS POINT
%     strain=B*U(sctrB);
%     stress=C*strain;
%     disp(e,1) = N'*U(2*sctr-1);
%     disp(e,2) = N'*U(2*sctr-0);
%     sigma(e,1)    = stress(1);
%     sigma(e,2)    = stress(3);
%     sigmaRef(e,1) = 1000/I0*(L-yPt(1))*yPt(2);
%     sigmaRef(e,2) = -1000/2/I0*(t^2/4-yPt(2)^2);
%     xcoord(e,:)     = yPt;
% end   % of element loop
% 
% fem=load('fem.mat');

% colordef white
% figure,set (gcf,'Color','w')
% set(gca,'FontSize',14)
% hold on
% plot(xcoord(:,2),sigma(:,2),'-o','MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',10.5);
% plot(fem.xcoord(:,2),fem.sigma(:,2),'s','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',6.5);
% h=legend('sigmaxx-exact','sigmaxx-coupling','sigmaxy-exact','sigmaxy-coupling');
% xlabel('x')
% ylabel('stresses along the vertical line x=11.85')
% grid on
% 
% colordef white
% figure,set (gcf,'Color','w')
% set(gca,'FontSize',14)
% hold on
% plot(xcoord(:,2),disp(:,2),'-o','MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',10.5);
% plot(fem.xcoord(:,2),fem.disp(:,2),'s','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',6.5);
% h=legend('sigmaxx-exact','sigmaxx-coupling','sigmaxy-exact','sigmaxy-coupling');
% xlabel('x')
% ylabel('stresses along the vertical line x=11.85')
% grid on

%%

aa=numx1;
elems = aa:numx1:aa+numx1*(numy1-1);

disp       = zeros(length(elems),2);
sigma      = zeros(length(elems),2);
sigmaRef   = zeros(length(elems),2);
xcoord     = zeros(length(elems),2);

for e=1:length(elems)
    ie = elems(e);
    sctr=element1(ie,:);
    nn   = length(sctr);
    sctrB(1:2:2*nn) = 2*sctr-1;
    sctrB(2:2:2*nn) = 2*sctr-0;
    nn=length(sctr);
    
    pt=[0 0];
    [N,dNdxi]=lagrange_basis(elemType,pt);
    J0=node1(sctr,:)'*dNdxi;
    invJ0=inv(J0);
    dNdx=dNdxi*invJ0;
    yPt=N'*node1(sctr,:);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % COMPUTE B MATRIX
    B(1,1:2:2*nn)  = dNdx(:,1)';
    B(2,2:2:2*nn)  = dNdx(:,2)';
    B(3,1:2:2*nn)  = dNdx(:,2)';
    B(3,2:2:2*nn)  = dNdx(:,1)';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % COMPUTE ELEMENT STRAIN AND STRESS AT STRESS POINT
    strain=B*U(sctrB);
    stress=C*strain;
    disp(e,1) = N'*U(2*sctr-1);
    disp(e,2) = N'*U(2*sctr-0);
    sigma(e,1)    = stress(1);
    sigma(e,2)    = stress(3);
    sigmaRef(e,1) = 1000/I0*(L-yPt(1))*yPt(2);
    sigmaRef(e,2) = -1000/2/I0*(t^2/4-yPt(2)^2);
    xcoord(e,:)     = yPt;
end   % of element loop

% fem=load('fem1.mat');
% 
% colordef white
% figure,set (gcf,'Color','w')
% set(gca,'FontSize',14)
% hold on
% plot(xcoord(:,2),sigma(:,2),'-o','MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',10.5);
% plot(fem.xcoord(:,2),sigmaRef(:,2),'s','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',6.5);
% h=legend('sigmaxx-exact','sigmaxx-coupling','sigmaxy-exact','sigmaxy-coupling');
% xlabel('x')
% ylabel('stresses along the vertical line x=11.85')
% grid on

colordef white
figure,set (gcf,'Color','w')
set(gca,'FontSize',14)
hold on
plot(xcoord(:,2),sigmaRef(:,2),'-o','MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',10.5);
plot(xcoord(:,2),sigma(:,2),'s','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',6.5);
h=legend('sigmaxx-exact','sigmaxx-coupling','sigmaxy-exact','sigmaxy-coupling');
xlabel('x')
ylabel('stresses along the vertical line x=11.85')
grid on

% xx= fem.sigma(:,2) - sigma(:,2);
% xx= xx./fem.sigma(:,2);
% 
% xx2 = xx.^2;
% err = sqrt(sum(xx2));

%% plot stresses at the top row elements
% GP=center point of each element

% elems = numx1*(numy1-1)+1:numx1*numy1;
% 
% sigma      = zeros(length(elems),2);
% sigmaRef   = zeros(length(elems),2);
% xcoord     = zeros(length(elems),2);
% 
% for e=1:length(elems)
%     ie = elems(e);
%     sctr=element1(ie,:);
%     nn   = length(sctr);
%     sctrB(1:2:2*nn) = 2*sctr-1;
%     sctrB(2:2:2*nn) = 2*sctr-0;
%     nn=length(sctr);
%     
%     pt=[0 0];
%     [N,dNdxi]=lagrange_basis(elemType,pt);
%     J0=node1(sctr,:)'*dNdxi;
%     invJ0=inv(J0);
%     dNdx=dNdxi*invJ0;
%     yPt=N'*node1(sctr,:);
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % COMPUTE B MATRIX
%     B(1,1:2:2*nn)  = dNdx(:,1)';
%     B(2,2:2:2*nn)  = dNdx(:,2)';
%     B(3,1:2:2*nn)  = dNdx(:,2)';
%     B(3,2:2:2*nn)  = dNdx(:,1)';
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % COMPUTE ELEMENT STRAIN AND STRESS AT STRESS POINT
%     strain=B*U(sctrB);
%     stress=C*strain;
%     sigma(e,1)    = stress(1);
%     sigma(e,2)    = stress(3);
%     sigmaRef(e,1) = 1000/I0*(L-yPt(1))*yPt(2);
%     sigmaRef(e,2) = -1000/2/I0*(t^2/4-yPt(2)^2);
%     xcoord(e,:)     = yPt;
% end   % of element loop
% 
% colordef white
% figure,set (gcf,'Color','w')
% set(gca,'FontSize',14)
% hold on
% plot(xcoord(:,1),sigmaRef(:,1),'k-','LineWidth',1.4);
% plot(xcoord(:,1),sigma(:,1),'o','MarkerEdgeColor','k',...
%     'MarkerFaceColor','g',...
%     'MarkerSize',6.5);
% plot(xcoord(:,1),sigmaRef(:,2),'k-','LineWidth',1.4);
% plot(xcoord(:,1),sigma(:,2),'s','MarkerEdgeColor','k',...
%     'MarkerFaceColor','g',...
%     'MarkerSize',6.5);
% h=legend('sigmaxx-exact','sigmaxx-coupling','sigmaxy-exact','sigmaxy-coupling');
% xlabel('x')
% ylabel('stresses along y=5.4')
% grid on
% %axis([0 5 -0.55 0])
% 
% %% plot stresses at the midline of the beam
% 
% no    = floor(numy1/2);
% elems = numx1*(no)+1:numx1*(no+1);
% 
% sigma      = zeros(length(elems),2);
% sigmaRef   = zeros(length(elems),2);
% xcoord     = zeros(length(elems),2);
% 
% for e=1:length(elems)
%     ie = elems(e);
%     sctr=element1(ie,:);
%     nn   = length(sctr);
%     sctrB(1:2:2*nn) = 2*sctr-1;
%     sctrB(2:2:2*nn) = 2*sctr-0;
%     nn=length(sctr);
%     
%     pt=[0 0];
%     [N,dNdxi]=lagrange_basis(elemType,pt);
%     J0=node1(sctr,:)'*dNdxi;
%     invJ0=inv(J0);
%     dNdx=dNdxi*invJ0;
%     yPt=N'*node1(sctr,:);
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % COMPUTE B MATRIX
%     B(1,1:2:2*nn)  = dNdx(:,1)';
%     B(2,2:2:2*nn)  = dNdx(:,2)';
%     B(3,1:2:2*nn)  = dNdx(:,2)';
%     B(3,2:2:2*nn)  = dNdx(:,1)';
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % COMPUTE ELEMENT STRAIN AND STRESS AT STRESS POINT
%     strain=B*U(sctrB);
%     stress=C*strain;
%     sigma(e,1)    = stress(1);
%     sigma(e,2)    = stress(3);
%     sigmaRef(e,1) = 1000/I0*(L-yPt(1))*yPt(2);
%     sigmaRef(e,2) = -1000/2/I0*(t^2/4-yPt(2)^2);
%     xcoord(e,:)     = yPt;
% end   % of element loop
% 
% colordef white
% figure,set (gcf,'Color','w')
% set(gca,'FontSize',14)
% hold on
% plot(xcoord(:,1),sigmaRef(:,1),'k-','LineWidth',1.4);
% plot(xcoord(:,1),sigma(:,1),'o','MarkerEdgeColor','k',...
%     'MarkerFaceColor','g',...
%     'MarkerSize',6.5);
% plot(xcoord(:,1),sigmaRef(:,2),'k--','LineWidth',1.4);
% plot(xcoord(:,1),sigma(:,2),'s','MarkerEdgeColor','k',...
%     'MarkerFaceColor','g',...
%     'MarkerSize',6.5);
% h=legend('sigmaxx-exact','sigmaxx-coupling','sigmaxy-exact','sigmaxy-coupling');
% xlabel('x')
% ylabel('stresses along y=0.3')
% grid on
% %axis([0 5 -0.55 0])

%% plot stresses at the coupling interface

% [W1,Q1]=quadrature( 3, 'GAUSS', 1 ); % two point quadrature
% 
% yy     = [];
% sXY    = [];
% sXYRef = [];
% 
% for e=1:numy1
%     sctrEdge = bndEdge1(e,:);
%     sctr1    = element1(bndMesh1(e),:);
%     pts1     = node1(sctr1,:);
%     nn   = length(sctr1);
%     sctrB(1:2:2*nn) = 2*sctr1-1;
%     sctrB(2:2:2*nn) = 2*sctr1-0;
%     for q=1:size(W1)
%         pt=Q1(q,:);
%         wt=W1(q);
%         [N,dNdxi]=lagrange_basis('L2',pt);  % element shape functions
%         J0=dNdxi'*node1(sctrEdge,:);
%         detJ0=norm(J0);
%         J0 = J0/detJ0;
%         
%         x=N'*node1(sctrEdge,:);
%         X1 = global2LocalMap(x,pts1,elemType);
%         
%         [N,dNdxi]=lagrange_basis(elemType,X1);
%         J0=node1(sctr,:)'*dNdxi;
%         invJ0=inv(J0);
%         dNdx=dNdxi*invJ0;
%         
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % COMPUTE B MATRIX
%         B(1,1:2:2*nn)  = dNdx(:,1)';
%         B(2,2:2:2*nn)  = dNdx(:,2)';
%         B(3,1:2:2*nn)  = dNdx(:,2)';
%         B(3,2:2:2*nn)  = dNdx(:,1)';
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % COMPUTE ELEMENT STRAIN AND STRESS AT STRESS POINT
%         strain=B*U(sctrB);
%         stress=C*strain;
%         yy = [yy;x(2)];
%         sXY = [sXY;stress(3)];
%         sXYRef = [sXYRef;-1000/2/I0*(t^2/4-x(2)^2)];
%     end
% end
% 
% colordef white
% figure,set (gcf,'Color','w')
% set(gca,'FontSize',14)
% hold on
% plot(yy,sXYRef,'k-','LineWidth',1.4);
% plot(yy,sXY,'o','MarkerEdgeColor','k',...
%     'MarkerFaceColor','g',...
%     'MarkerSize',6.5);
% h=legend('exact','coupling');
% xlabel('y')
% ylabel('stresses along the coupling interface')
% grid on

%% contour plot of stress field

stress=zeros(size(element1,1),size(element1,2),3);

stressPoints=[-1 -1;1 -1;1 1;-1 1];


for e=1:size(element1,1)
    sctr=element1(e,:);
    sctrB(1:2:2*nn) = 2*sctr-1;
    sctrB(2:2:2*nn) = 2*sctr-0;
    nn=length(sctr);
    Ce=C;
    
    for q=1:nn
        pt=stressPoints(q,:);
        [N,dNdxi]=lagrange_basis(elemType,pt);
        J0=node1(sctr,:)'*dNdxi;
        invJ0=inv(J0);
        dNdx=dNdxi*invJ0;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % COMPUTE B MATRIX
        B(1,1:2:2*nn)  = dNdx(:,1)';
        B(2,2:2:2*nn)  = dNdx(:,2)';
        B(3,1:2:2*nn)  = dNdx(:,2)';
        B(3,2:2:2*nn)  = dNdx(:,1)';
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % COMPUTE ELEMENT STRAIN AND STRESS AT STRESS POINT
        strain=B*U(sctrB);
        sigma =C*strain;
        stress(e,q,1:3)=sigma;
        stress(e,q,4)  = sqrt(sigma(1)^2+sigma(2)^2-sigma(1)*sigma(2)+3*(sigma(3)^2));
    end
end   % of element loop


%print(gcf, '-depsc', '-painters', 'output.eps')
%print(fn,?-djpeg90?,[?beam_?,elemType,?_sigma?,num2str(stressComp),?.jpg?])

numNode  = size(node1,1);

% normal stresses
sigmaXX = zeros(numNode,2);
sigmaYY = zeros(numNode,2);
sigmaXY = zeros(numNode,2);
sigmaVM = zeros(numNode,2);

for e=1:size(element1,1)
    connect = element1(e,:);
    for in=1:4
        nid = connect(in);
        sigmaXX(nid,:) = sigmaXX(nid,:) + [stress(e,in,1) 1];
        sigmaYY(nid,:) = sigmaYY(nid,:) + [stress(e,in,2) 1];
        sigmaXY(nid,:) = sigmaXY(nid,:) + [stress(e,in,3) 1];
        sigmaVM(nid,:) = sigmaVM(nid,:) + [stress(e,in,4) 1];
    end
end

% Average nodal stress values (learned from Mathiew Pais XFEM code)
sigmaXX(:,1) = sigmaXX(:,1)./sigmaXX(:,2); sigmaXX(:,2) = [];
sigmaYY(:,1) = sigmaYY(:,1)./sigmaYY(:,2); sigmaYY(:,2) = [];
sigmaXY(:,1) = sigmaXY(:,1)./sigmaXY(:,2); sigmaXY(:,2) = [];
sigmaVM(:,1) = sigmaVM(:,1)./sigmaVM(:,2); sigmaVM(:,2) = [];
%
% plotStress(sigmaXX,sigmaXY,sigmaYY,element1,node1);

scaleFact=100;

figure
clf
plot_field(node1+scaleFact*[Ux1 Uy1],element1,elemType,sigmaVM);
plot_mesh(node2+scaleFact*[zeros(numnode2,1) Ux2],element2,elemTypeB,'b-',2);
%plot_mesh(node2+scaleFact*[Ux2 Uy2],element2,elemType,'r.-',1);
hold on
%plot_mesh(node+scaleFact*[U(xs) U(ys)],element,elemType,'g.-');
%plot_mesh(node,element,elemType,'w--');
colorbar
axis off
%%
