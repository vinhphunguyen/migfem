% This file implements the Arlequin method.
% With coarse coupling and H1 norm and constant weight functions.
% Only compatible and structured mesh is supported.
%
% Used to check the implementation of the Arlequin method.
%
% Vinh Phu Nguyen
% Ton Duc Thang University, HCMC, Vietnam
% November 2012

addpath ../fem_util/
addpath ../gmshFiles/
addpath ../post-processing/
addpath ../fem-functions/
addpath ../analytical-solutions/

clear all
colordef black
state = 0;
tic;

global elemType numdofs C alpha1 alpha2 ne lagDom l couDom2 fine2VirtualMap

% MATERIAL PROPERTIES
E0  = 30e6;  % Young?s modulus
nu0 = 0.3;  % Poisson?s ratio

% BEAM PROPERTIES
W  = 10;     % 
L  = 20;

plotMesh  = 1;

% LOAD

sigmato = 100;

% COMPUTE ELASTICITY MATRIX
% Plane stress

C=E0/(1-nu0^2)*[   1      nu0          0;
    nu0        1          0;
    0        0  (1-nu0)/2 ];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE FINITE ELEMENT MESH
%

disp([num2str(toc),'   GENERATING MESH'])

% MESH PROPERTIES

elemType = 'Q4'; % the element type used in the FEM simulation;

% domain 1
numy1     = 16;
numx1     = 8;

% size of coupling zone

% meshing for domain1 (fine mesh)
nnx=numx1+1;
nny=numy1+1;
node=square_node_array([0 0],[W 0],[W L],[0 L],nnx,nny);
inc_u=1;
inc_v=nnx;
node_pattern=[ 1 2 nnx+2 nnx+1 ];
element =make_elem(node_pattern,numx1,numy1,inc_u,inc_v);


% PLOT MESH
if ( plotMesh )
    clf
    hold on
    %plot_mesh(node1b,element1b,elemType,'g.-',3.5);
    %plot_mesh(node1a,element1a,elemType,'g.-',3.5);    
    plot_mesh(node,element,elemType,'g.-',3.5);              
    axis off
    %axis([0 L -c c])
end

% DEFINE BOUNDARIES

edgeElemType='L2';

% GET NODES ON DISPLACEMENT BOUNDARY
%      Here we get the nodes on the essential boundaries

fixedNode=find(node(:,2)==0);
topNode  =find(node(:,2)==L);

topEdge   = zeros(numx1,2);

for i=1:numx1
    topEdge(i,:) = topNode(i:i+1);
end

uFixed=zeros(1,length(fixedNode))';  % a vector of u_x for the nodes
vFixed=zeros(1,length(fixedNode))';

numdofs = size(node,1);

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

[W,Q]=quadrature( 2, 'GAUSS', 2 ); % 2x2 Gaussian quadrature

for e=1:size(element,1)                          % start of element loop    
    [Ke,sctrB]     = stiffnessMatrix (e,element,node,W,Q);
    K(sctrB,sctrB) = K(sctrB,sctrB) + Ke;
end


% Compute external force vector as usual

[W,Q]=quadrature( 3, 'GAUSS', 1 ); % three point quadrature

% TOP EDGE
for e=1:size(topEdge,1) % loop over the elements in the right edge
    sctr=topEdge(e,:);  % scatter vector for the element
    sctrx=sctr;           % x scatter vector
    sctry=sctrx+numdofs;  % y scatter vector
    for q=1:size(W,1)
        pt=Q(q,:);
        wt=W(q);
        [N,dNdxi]=lagrange_basis(edgeElemType,pt);  % element shape functions
        J0=dNdxi'*node(sctr,:);
        detJ0=norm(J0);                
        f(sctrx)=f(sctrx)+N*sigmato*detJ0*wt;
    end % of quadrature loop
end  % of element loop


%%%%%%%%%%%%%%%%%%% END OF STIFFNESS MATRIX COMPUTATION %%%%%%%%%%%%%%%%%%%


% APPLY ESSENTIAL BOUNDARY CONDITIONS
disp([num2str(toc),'   APPLYING BOUNDARY CONDITIONS'])
bcwt=mean(diag(K)); % a measure of the average size of an element in K
% used to keep the conditioning of the K matrix
udofs=fixedNode;           % global indecies of the fixed x displacements
vdofs=fixedNode+numdofs;   % global indecies of the fixed y displacements
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

Ux = U(1:size(node,1));
Uy = U([1:size(node,1)]+numdofs);

% Here we plot the stresses and displacements of the solution. As with the
% mesh generation section we don?t go into too much detail - use help
% ?function name? to get more details.
disp([num2str(toc),'   POST-PROCESSING'])

scaleFact=1e3;
fn=2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT DEFORMED DISPLACEMENT PLOT
%Ux = U(xs);
%Uy = U(ys);

figure(fn)
clf
%plot_field(node1+scaleFact*[Ux Uy],element1,elemType,U(ys));
%plot_field(node2+scaleFact*[Ux Uy],element2,elemType,U(ys));

hold on
plot_field(node+scaleFact*[Ux Uy],element,elemType,Uy);
plot_mesh(node+scaleFact*[Ux Uy],element,elemType,'r.-',1.1);
colorbar
fn=fn+1;
title('DEFORMED DISPLACEMENT IN Y-DIRECTION')

  
% COMPUTE STRESS

stress=zeros(size(element1,1),size(element1,2),3);

stressPoints=[-1 -1;1 -1;1 1;-1 1];


for e=1:size(element1,1)
    sctr=element1(e,:);
    sctrB=[sctr sctr+numdofs];
    nn=length(sctr);
    for q=1:nn
        pt=stressPoints(q,:);
        [N,dNdxi]=lagrange_basis(elemType,pt);
        J0=node1(sctr,:)'*dNdxi;
        invJ0=inv(J0);
        dNdx=dNdxi*invJ0;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % COMPUTE B MATRIX
        B=zeros(3,2*nn);
        B(1,1:nn)       = dNdx(:,1)';
        B(2,nn+1:2*nn)  = dNdx(:,2)';
        B(3,1:nn)       = dNdx(:,2)';
        B(3,nn+1:2*nn)  = dNdx(:,1)';
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % COMPUTE ELEMENT STRAIN AND STRESS AT STRESS POINT
        strain=B*U(sctrB);
        stress(e,q,:)=C*strain;
    end
end   % of element loop

stressComp=1;
figure(fn)
clf
plot_field(node1+scaleFact*[Ux1 Uy1],element1,elemType,stress(:,:,stressComp));
plot_mesh(node2+scaleFact*[Ux2 Uy2],element2,elemType,'r.-');
hold on
%plot_mesh(node+scaleFact*[U(xs) U(ys)],element,elemType,'g.-');
%plot_mesh(node,element,elemType,'w--');
colorbar
fn=fn+1;
title('DEFORMED STRESS PLOT, BENDING COMPONENT')
