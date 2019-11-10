% This file implements the Nitsche method to join two mechanical models:
% a 2D continuum model and a frame model.
% The beam elements can be placed arbitrarily in space.
% This is a prototype for solving frame problems with continuum elements
% in region of high gradients and beam elements in other places.
%
% Dofs: 
%     continuum nodes: ux and uy
%     beam nodes:      ux,uy and rotation
% Assembly: assume that continuum nodes are stored first followed by beam nodes.
%
% Discretisation: standard Lagrange finite elements.
%
% Problem: Timoshenko beam in bending.
% 
% Vinh Phu Nguyen
% Cardiff University, UK
% 7 June 2013

addpath ../fem_util/
addpath ../gmshFiles/
addpath ../post-processing/
addpath ../fem-functions/
addpath ../analytical-solutions/


clear all
colordef white
state = 0;
tic;

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,'splinecurve.eps',opts)


% MATERIAL PROPERTIES
E0  = 1000;  % Young?s modulus
nu0 = 0.3;  % Poisson?s ratio

% BEAM PROPERTIES
L  = 8;        % length of the beam
c  = 0.5;      % the distance of the outer fiber of the beam from the mid-line
t  = 2*c;      % thickness, area = bxt
b  = 1; 

I0 = 2*c^3/3;  % the second polar moment of inertia of the beam cross-section.
q  = 0;        % uniformly distributed loads
F  = 1;        % end force
G  = E0/2/(1+nu0);
k  = E0*I0; 
A  = b*t;
shearFactor = 5/6;

plotMesh  = 1;

% TIP LOAD
P = 1000; % the peak magnitude of the traction at the right edge

stressState = 'PLANE_STRESS';
% COMPUTE ELASTICITY MATRIX
C = elasticityMatrix(E0,nu0,stressState);

noGPs = 2;

% penalty parameter

alpha = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE FINITE ELEMENT MESH
%

disp([num2str(toc),'   GENERATING MESH'])

% MESH PROPERTIES

elemType  = 'Q4'; % the element type used in the FEM simulation;
elemTypeB = 'L2'; % the element type used in the FEM simulation;

% domain 1
numx1     = 14;
numy1     = 28;

L1 = 4.;
% meshing for domain1
nnx=numx1+1;
nny=numy1+1;
%node1=square_node_array([0 0],[t 0],[t L1],[0 L1],nnx,nny);
node1=square_node_array([-c 0],[c 0],[c L],[-c L],nnx,nny);
%node1=square_node_array([-t 0],[0 0],[0 L1],[-t L1],nnx,nny);
inc_u=1;
inc_v=nnx;
node_pattern=[ 1 2 nnx+2 nnx+1 ];
element1 =make_elem(node_pattern,numx1,numy1,inc_u,inc_v);



% PLOT MESH
if ( plotMesh )  
    clf
    plot_mesh(node1,element1,elemType,'g.-',2.4);    
    hold on
    axis off
    %axis([0 L -c c])
end



% DEFINE BOUNDARIES

edgeElemType='L2';

% GET NODES ON DISPLACEMENT BOUNDARY
%      Here we get the nodes on the essential boundaries

fixedNode = find(node1(:,2)==0);  
midNode1  = find(node1(:,1)==0);  
rightNode = find(abs(node1(:,2)-L)<1e-14); 

uFixed=zeros(1,length(fixedNode))';  % a vector of u_x for the nodes
vFixed=zeros(1,length(fixedNode))';     

topEdge  = zeros(numx1,2);

for i=1:numx1
    topEdge(i,:) = rightNode(i:i+1);
end

numdofs = 2*size(node1,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE SYSTEM DATA STRUCTURES

disp([num2str(toc),' INITIALIZING DATA STRUCTURES'])

f=zeros(numdofs,1);          % external load vector
K=zeros(numdofs,numdofs); % stiffness matrix

% ******************************************************************************
% ***                          P R O C E S S I N G                           ***
% ******************************************************************************
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%% COMPUTE STIFFNESS MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp([num2str(toc),'   COMPUTING STIFFNESS MATRIX'])

[W,Q]=quadrature( 2, 'GAUSS', 2 ); % 2x2 Gaussian quadrature

% for domain 1

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
        % B matrix
        B(1,1:2:2*nn)  = dNdx(:,1)';
        B(2,2:2:2*nn)  = dNdx(:,2)';
        B(3,1:2:2*nn)  = dNdx(:,2)';
        B(3,2:2:2*nn)  = dNdx(:,1)';        
        % COMPUTE ELEMENT STIFFNESS AT QUADRATURE POINT
        K(sctrB,sctrB)=K(sctrB,sctrB)+B'*C*B*W(q)*det(J0);
    end  % of quadrature loop
end    

[W,Q]=quadrature( 2, 'GAUSS', 1 ); % 2x2 Gaussian quadrature

% RIGHT EDGE
for e=1:size(topEdge,1) % loop over the elements in the right edge
    sctr=topEdge(e,:);  % scatter vector for the element
    sctrx=2*sctr-1;           % x scatter vector
    
    for q=1:size(W,1)
        pt=Q(q,:);
        wt=W(q);
        [N,dNdxi]=lagrange_basis(edgeElemType,pt);  % element shape functions
        J0=dNdxi'*node1(sctr,:);
        detJ0=norm(J0);                                
        fyPt=F;
        f(sctrx)=f(sctrx)+N*fyPt*detJ0*wt;  %tor        
    end % of quadrature loop
end  % of element loop



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
f(vdofs)=bcwt*speye(length(udofs))*vFixed;

% SOLVE SYSTEM
disp([num2str(toc),'   SOLVING SYSTEM'])
U=K\f;

%******************************************************************************
%*** POST - PROCESSING *** 
%***************************************************

numnode1 = size(node1,1);
index1 = 1:numnode1;
Ux1  = U(2*index1-1);
Uy1  = U(2*index1);


% Here we plot the stresses and displacements of the solution. As with the
% mesh generation section we don?t go into too much detail - use help
% ?function name? to get more details.
disp([num2str(toc),'   POST-PROCESSING'])

scaleFact=3.;
fn=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT DEFORMED DISPLACEMENT PLOT

colordef black
figure
clf
hold on
%plot_field(node1+scaleFact*[Ux1 Uy1],element1,elemType,Ux1);
%plot_field(node2+scaleFact*[Ux2 Uy2],element2,elemType,Ux2);
plot_mesh(node1+scaleFact*[Ux1 Uy1],element1,elemType,'w.-',1);
%colorbar
axis off
%title('DEFORMED DISPLACEMENT IN Y-DIRECTION')

% Comapre numerical displacement to exact value
 
Uxm1  = U(2*midNode1-1);
Uym1  = U(2*midNode1);

xx = [node1(midNode1,2)];
u  = [-Uxm1];


% exact solution u = 

x      = linspace(0,L,400);
uExact = -F/(6*k)*(3*L*x.^2-x.^3);

% solution from 2D continuum elements

uCont     = csvread('beam-2D.csv');
uCoup     = csvread('beam-1D.csv');

colordef white
figure,set (gcf,'Color','w')
hold on
plot(x,uExact,'r-','LineWidth',1.4);
plot(xx,u,'k-','MarkerSize',12,'LineWidth',1.4);
plot(uCont(:,1),uCont(:,2),'b--','LineWidth',1.4);
plot(uCoup(:,1),uCoup(:,2),'cy.-','LineWidth',1.4);
h=legend('exact','coupling','continuum','beam');
xlabel('x')
ylabel('w')
%axis([0 5 -0.55 0])

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
        stress(e,q,:)=Ce*strain;
    end
end   % of element loop


% sigmaXX = zeros(size(node1,1),2);
% sigmaYY = zeros(size(node1,1),2);
% sigmaXY = zeros(size(node1,1),2);
% 
% for e=1:size(element1,1)
%     connect = element1(e,:);
%     for in=1:4
%         nid = connect(in);
%         sigmaXX(nid,:) = sigmaXX(nid,:) + [stress(e,in,1) 1];
%         sigmaYY(nid,:) = sigmaYY(nid,:) + [stress(e,in,2) 1];
%         sigmaXY(nid,:) = sigmaXY(nid,:) + [stress(e,in,3) 1];               
%     end
% end
% 
% % Average nodal stress values (learned from Mathiew Pais XFEM code)
% sigmaXX(:,1) = sigmaXX(:,1)./sigmaXX(:,2); sigmaXX(:,2) = [];
% sigmaYY(:,1) = sigmaYY(:,1)./sigmaYY(:,2); sigmaYY(:,2) = [];
% sigmaXY(:,1) = sigmaXY(:,1)./sigmaXY(:,2); sigmaXY(:,2) = [];
% 
% plotStress(sigmaXX,sigmaXY,sigmaYY,element1,node1);

stressComp=2;
figure
clf
plot_field(node1+scaleFact*[Ux1 Uy1],element1,elemType,stress(:,:,stressComp));
%plot_mesh(node2+scaleFact*[Ux2 Uy2],element2,elemType,'r.-',1);
hold on
%plot_mesh(node+scaleFact*[U(xs) U(ys)],element,elemType,'g.-');
%plot_mesh(node,element,elemType,'w--');
colorbar

