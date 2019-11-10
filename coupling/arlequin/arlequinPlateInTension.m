% This file implements the Arlequin method.
% With coarse coupling and H1 norm and constant weight functions.
% Only compatible and structured mesh is supported.
%
% Plate in tension. Two ends of the plate are meshed with coarse mesh
% and the middle part is discretize with a refined mesh. There are two
% coupling zones.
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
H1 = 5;
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
numy1     = 8;
numx1     = 8;

% size of coupling zone

L1 = 8;
dx = 6*(L1/numy1);
l  = dx/6; l=2;

% domain 2
numy2     = 12;
numx2     = 16;

% meshing for domain1 (fine mesh)
nnx=numx2+1;
nny=numy2+1;
node1=square_node_array([0 7],[W 7],[W 13],[0 13],nnx,nny);
inc_u=1;
inc_v=nnx;
node_pattern=[ 1 2 nnx+2 nnx+1 ];
element1 =make_elem(node_pattern,numx2,numy2,inc_u,inc_v);

% meshing for domain2 (coarse mesh)
nnx=numx1+1;
nny=numy1+1;
node2a=square_node_array([0 0],[W 0],[W 8],[0 8],nnx,nny);
inc_u=1;
inc_v=nnx;
node_pattern=[ 1 2 nnx+2 nnx+1 ];
element2a =make_elem(node_pattern,numx1,numy1,inc_u,inc_v);

node2b=square_node_array([0 12],[W 12],[W 20],[0 20],nnx,nny);
inc_u=1;
inc_v=nnx;
node_pattern=[ 1 2 nnx+2 nnx+1 ];
element2b =make_elem(node_pattern,numx1,numy1,inc_u,inc_v);

% combine these two meshes to make one mesh

node2=[node2a;node2b];
element2=[element2a;element2b+size(node2a,1)];


couDom1 = [];
couDom2 = [];

matDom1 = ones(size(element1,1),1);
matDom2 = ones(size(element2,1),1);

keys = [];
vals = [];

for e=1:size(element1,1)
    sctr   = element1(e,:);
    nodes  = node1(sctr,:);
    center = mean(nodes);
    
    if center(2) < 8 || center(2) > 12
        couDom1 = [couDom1 e];
        matDom2(e) = 2;
        keys = [keys e];
    end
end

for e=1:size(element2a,1)
    sctr   = element2a(e,:);
    nodes  = node2a(sctr,:);
    center = mean(nodes);
    
    if center(2) > 8-2
        couDom2 = [couDom2 e];
        matDom1(e) = 2;
        %keys = [keys e];
        %vals = [vals 1];
    end
end

for e=1:size(element2b,1)
    sctr   = element2b(e,:);
    nodes  = node2b(sctr,:);
    center = mean(nodes);
    
    if center(2) < 14
        couDom2 = [couDom2 e+size(element2a,1)];
        matDom1(e) = 2;
        %keys = [keys e];
        %vals = [vals 1];
    end
end

% THIS IS PROBLEM DEPENDENT !!!

vals=[1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8 ...
      1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8 ...
      9 9 10 10 11 11 12 12 13 13 14 14 15 15 16 16 ...
      9 9 10 10 11 11 12 12 13 13 14 14 15 15 16 16];% ...
%       17 17 18 18 19 19 20 20 21 21 22 22 23 23 24 24 ...
%       17 17 18 18 19 19 20 20 21 21 22 22 23 23 24 24 ...
%       25 25 26 26 27 27 28 28 29 29 30 30 31 31 32 32 ...
%       25 25 26 26 27 27 28 28 29 29 30 30 31 31 32 32 ...
%       ];



fine2VirtualMap = containers.Map(keys,vals);

% virtual elements for discretizing Lagrange multipliers
% coarse couping => mesh for Lagrange multipliers coincide
% the coarse mesh in the coupling zone

%lagDom  = [1 2 4 3; 3 4 6 5];

numxL = numx1 ;
numyL = 4;

nnx=numxL+1;
nny=numyL+1;
inc_u=1;
inc_v=nnx;
node_pattern=[ 1 2 nnx+2 nnx+1 ];
lagDom =make_elem(node_pattern,numxL,numyL,inc_u,inc_v);

lagDom  = lagDom + size(node1,1) +  size(node2,1);

% partition of unity for energy

alpha1 = 0.9;
alpha2 = 0.1;

% PLOT MESH
if ( plotMesh )
    clf
    hold on
    %plot_mesh(node1b,element1b,elemType,'g.-',3.5);
    %plot_mesh(node1a,element1a,elemType,'g.-',3.5);    
    plot_mesh(node1,element1,elemType,'g.-',1.11);  
    %plot(node1(:,1),node1(:,2),'ro')
    plot_mesh(node2,element2,elemType,'r-',3);    
    axis on
    %axis([0 L -c c])
end

% DEFINE BOUNDARIES

edgeElemType='L2';

% GET NODES ON DISPLACEMENT BOUNDARY
%      Here we get the nodes on the essential boundaries

fixedNode=find(node2(:,2)==0)+size(node1,1);
topNode  =find(node2b(:,2)==L);

topEdge   = zeros(numx1,2);

for i=1:numx1
    topEdge(i,:) = topNode(i:i+1);
end

uFixed=zeros(1,length(fixedNode))';  % a vector of u_x for the nodes
vFixed=zeros(1,length(fixedNode))';

numdofs = size(node1,1) + size(node2,1) + ...
          length(unique(lagDom));

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

% for domain 1 (fine mesh)

for e=1:size(element1,1)                          % start of element loop
    % skip elements in the coupling domain
    if ~ismember(e,couDom1)
        [Ke,sctrB]     = stiffnessMatrix (e,element1,node1,W,Q);
        K(sctrB,sctrB) = K(sctrB,sctrB) + Ke;
    else
        [Ke,Ce,sctrB,sctrLB] = stiffnessCoupledMatrix...
            (e,element1,element2,node1,node2,W,Q);
        
        K(sctrB,sctrB)  = K(sctrB,sctrB)  + Ke;
        K(sctrB,sctrLB) = K(sctrB,sctrLB) - Ce';
        K(sctrLB,sctrB) = K(sctrLB,sctrB) - Ce;
    end
end

% for domain 2

%element2 = element2 + size(node1,1);

for e=1:size(element2,1)                          
    e;
    if ~ismember(e,couDom2)
        sctr  = element2(e,:);
        sctr2 = sctr + size(node1,1);
                
        [Ke,sctrB] = stiffnessMatrix (e,element2,node2,W,Q);
        sctrB = [ sctr2 sctr2+numdofs ]; % vector that scatters a B matrix
        K(sctrB,sctrB) = K(sctrB,sctrB) + Ke;
    else
        sctr  = element2(e,:);
        sctr2 = sctr + size(node1,1);
        sctrL  = lagDom(find(couDom2==e),:);
        sctrB  = [ sctr2 sctr2+numdofs ];
        sctrLB = [ sctrL sctrL+numdofs ];
        nn=length(sctr);
        
        Ke = zeros(8,8);
        Ce = zeros(8,8);
        
        for q=1:size(W,1)
            pt=Q(q,:);
            wt=W(q);
            [N,dNdxi]=lagrange_basis(elemType,pt); % element shape functions
            J0=node2(sctr,:)'*dNdxi;                % element Jacobian matrix
            invJ0=inv(J0);
            dNdx=dNdxi*invJ0;
            
            B(1,1:nn)       = dNdx(:,1)';
            B(2,nn+1:2*nn)  = dNdx(:,2)';
            B(3,1:nn)       = dNdx(:,2)';
            B(3,nn+1:2*nn)  = dNdx(:,1)';
            
            Nm(1,1:nn)       = N';
            Nm(2,nn+1:2*nn)  = N';
            
            % COMPUTE ELEMENT STIFFNESS AT QUADRATURE POINT
            
            Ke = Ke + alpha2*B'*C*B*W(q)*det(J0);
            Ce = Ce + Nm'*Nm*W(q)*det(J0) + B'*B*W(q)*det(J0)*l*l;
        end  % of quadrature loop
        Ke;
        Ce;
        K(sctrB,sctrB)  = K(sctrB,sctrB)  + Ke;
        K(sctrB,sctrLB) = K(sctrB,sctrLB) + Ce';
        K(sctrLB,sctrB) = K(sctrLB,sctrB) + Ce;
    end
end

% Compute external force vector as usual

[W,Q]=quadrature( 2, 'GAUSS', 1 ); % three point quadrature

% TOP EDGE
for e=1:size(topEdge,1) % loop over the elements in the right edge
    sctr=topEdge(e,:);  % scatter vector for the element
    sctrx=sctr+ size(node2a,1) + size(node1,1);           % x scatter vector
    sctry=sctrx+numdofs;  % y scatter vector
    for q=1:size(W,1)
        pt=Q(q,:);
        wt=W(q);
        [N,dNdxi]=lagrange_basis(edgeElemType,pt);  % element shape functions
        J0=dNdxi'*node2b(sctr,:);
        detJ0=norm(J0);                
        f(sctrx)=f(sctrx)-N*sigmato*detJ0*wt;
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

Ux1 = U(1:size(node1,1));
Uy1 = U([1:size(node1,1)]+numdofs);

Ux2 = U(1+size(node1,1):size(node1,1)+size(node2,1));
Uy2 = U([1+size(node1,1):size(node1,1)+size(node2,1)]+numdofs);

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
plot_mesh(node1+scaleFact*[Ux1 Uy1],element1,elemType,'g.-',1.2);
plot_mesh(node2+scaleFact*[Ux2 Uy2],element2,elemType,'r.-',1.2);
plot_field(node1+scaleFact*[Ux1 Uy1],element1,elemType,Uy1);
plot_field(node2+scaleFact*[Ux2 Uy2],element2,elemType,Uy2);
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

% stressComp=1;
% figure(fn)
% clf
% plot_field(node1+scaleFact*[Ux1 Uy1],element1,elemType,stress(:,:,stressComp));
% plot_mesh(node2+scaleFact*[Ux2 Uy2],element2,elemType,'r.-',1.2);
% hold on
% %plot_mesh(node+scaleFact*[U(xs) U(ys)],element,elemType,'g.-');
% %plot_mesh(node,element,elemType,'w--');
% colorbar
% fn=fn+1;
% title('DEFORMED STRESS PLOT, BENDING COMPONENT')
