% This file implements the Arlequin method.
% With coarse coupling and H1 norm and constant weight functions.
% Only compatible and structured mesh is supported.
%
% Better implementation.
%
% Timoshenko beam in bending.
%
% Vinh Phu Nguyen
% Ton Duc Thang University, HCMC, Vietnam
% November 2012

addpath ../fem_util/
addpath ../gmshFiles/
addpath ../fem/
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
L  = 48;     % length of the beam
c  = 6;      % the distance of the outer fiber of the beam from the mid-line

plotMesh  = 1;

% TIP LOAD
P = 1000; % the peak magnitude of the traction at the right edge
I0=2*c^3/3;  % the second polar moment of inertia of the beam cross-section.


% COMPUTE ELASTICITY MATRIX
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
numx1     = 20;

% size of coupling zone

L1 = 20;
dx = 4*(L1/numx1);
l  = dx;

% domain 2
numy2     = 4;
numx2     = 15;
L =numx2*2+L1-dx;

% meshing for domain1
nnx=numx1+1;
nny=numy1+1;
node1=square_node_array([0 -c],[L1 -c],[L1 c],[0 c],nnx,nny);
inc_u=1;
inc_v=nnx;
node_pattern=[ 1 2 nnx+2 nnx+1 ];
element1 =make_elem(node_pattern,numx1,numy1,inc_u,inc_v);

% meshing for domain2
nnx=numx2+1;
nny=numy2+1;
node2=square_node_array([L1-dx -c],[L -c],[L c],[L1-dx c],nnx,nny);
inc_u=1;
inc_v=nnx;
node_pattern=[ 1 2 nnx+2 nnx+1 ];
element2 =make_elem(node_pattern,numx2,numy2,inc_u,inc_v);

no = 9;
%couDom1 = [no no+1 ...
%           no+numx1 no+numx1+1 ...
%           2*numx1+no 2*numx1+no+1 ...
%           3*numx1+no 3*numx1+no+1];

couDom1 = [];

matDom1 = ones(size(element1,1),1);
matDom2 = ones(size(element2,1),1);

keys = [];
vals = [];

for e=1:size(element1,1)
    sctr   = element1(e,:);
    nodes  = node1(sctr,:);
    center = mean(nodes);
    
    if center(1) > L1-dx
        couDom1 = [couDom1 e];
        matDom1(e) = 2;
        keys = [keys e];
        vals = [vals 1];
    end
end

% THIS IS PROBLEM DEPENDENT !!!

vals=[1 1 2 2 1 1 2 2 3 3 4 4 3 3 4 4 ...
    5 5 6 6 5 5 6 6 7 7 8 8 7 7 8 8];

fine2VirtualMap = containers.Map(keys,vals);

%couDom2 = [1 1+numx2];

couDom2 = [];

for e=1:size(element2,1)
    sctr   = element2(e,:);
    nodes  = node2(sctr,:);
    center = mean(nodes);
    
    if center(1) < L1
        couDom2 = [couDom2 e];
        matDom2(e) = 2;
    end
end

% virtual elements for discretizing Lagrange multipliers
% coarse couping => mesh for Lagrange multipliers coincide
% the coarse mesh in the coupling zone

%lagDom  = [1 2 4 3; 3 4 6 5];

numxL = 2 ;
numyL = 4;

nnx=numxL+1;
nny=numyL+1;
inc_u=1;
inc_v=nnx;
node_pattern=[ 1 2 nnx+2 nnx+1 ];
lagDom =make_elem(node_pattern,numxL,numyL,inc_u,inc_v);

lagDom  = lagDom + size(node1,1) + size(node2,1);

% partition of unity for energy

alpha1 = 0.7;
alpha2 = 0.3;

% PLOT MESH
if ( plotMesh )
    clf
    plot_mesh(node1,element1,elemType,'g.-',2.4);
    plot_mesh(node2,element2,elemType,'r.--',1);
    hold on
    axis off
    axis([0 L -c c])
end

% DEFINE BOUNDARIES

edgeElemType='L2';

% GET NODES ON DISPLACEMENT BOUNDARY
%      Here we get the nodes on the essential boundaries

fixedNode=find(node1(:,1)==0);
rightNode=find(node2(:,1)==L);

rightEdge   = zeros(numy2,2);

for i=1:numy2
    rightEdge(i,:) = rightNode(i:i+1);
end

uFixed=zeros(1,length(fixedNode))';  % a vector of u_x for the nodes
vFixed=zeros(1,length(fixedNode))';

numdofs = size(node1,1) + size(node2,1) + ...
    length(unique(element2(couDom2,:)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE SYSTEM DATA STRUCTURES

disp([num2str(toc),' INITIALIZING DATA STRUCTURES'])

f=zeros(2*numdofs,1);          % external load vector
K=zeros(2*numdofs,2*numdofs); % stiffness matrix


%xs=1:numnode;                  % x portion of u and v vectors
%ys=(numnode+1):2*numnode;      % y portion of u and v vectors

% ******************************************************************************
% ***                          P R O C E S S I N G                           ***
% ******************************************************************************
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%% COMPUTE STIFFNESS MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp([num2str(toc),'   COMPUTING STIFFNESS MATRIX'])

[W,Q]=quadrature( 2, 'GAUSS', 2 ); % 2x2 Gaussian quadrature

% for domain 1

i = 0;

for e=1:size(element1,1)                          % start of element loop
    % skip elements in the coupling domain
    if ~ismember(e,couDom1)
        [Ke,sctrB]     = stiffnessMatrix (e,element1,node1,W,Q);
        K(sctrB,sctrB) = K(sctrB,sctrB) + Ke;
    else
        i = i + 1;
        [Ke,Ce,sctrB,sctrLB] = stiffnessCoupledMatrix...
            (e,element1,element2,node1,node2,W,Q);
        
        K(sctrB,sctrB)  = K(sctrB,sctrB)  + Ke;
        K(sctrB,sctrLB) = K(sctrB,sctrLB) - Ce';
        K(sctrLB,sctrB) = K(sctrLB,sctrB) - Ce;
    end
end

% for domain 2

%element2 = element2 + size(node1,1);

for e=1:size(element2,1)                          % start of element loop
    if ~ismember(e,couDom2)
        sctr  = element2(e,:);           % element scatter vector
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
            
            BB = B;
            %BB(3,:) = 0.5 * BB(3,:);
            
            Nm(1,1:nn)       = N';
            Nm(2,nn+1:2*nn)  = N';
            
            % COMPUTE ELEMENT STIFFNESS AT QUADRATURE POINT
            
            Ke = Ke + alpha2*B'*C*B*W(q)*det(J0);
            Ce = Ce + Nm'*Nm*W(q)*det(J0) + B'*BB*W(q)*det(J0)*l*l;
        end  % of quadrature loop
        Ke;
        Ce;
        K(sctrB,sctrB)  = K(sctrB,sctrB)  + Ke;
        K(sctrB,sctrLB) = K(sctrB,sctrLB) + Ce';
        K(sctrLB,sctrB) = K(sctrLB,sctrB) + Ce;
    end
end

% Compute external force vector as usual

[W,Q]=quadrature( 3, 'GAUSS', 1 ); % three point quadrature

% RIGHT EDGE
for e=1:size(rightEdge,1) % loop over the elements in the right edge
    sctr=rightEdge(e,:);  % scatter vector for the element
    sctrx=sctr+ size(node1,1);           % x scatter vector
    sctry=sctrx+numdofs;  % y scatter vector
    for q=1:size(W,1)
        pt=Q(q,:);
        wt=W(q);
        [N,dNdxi]=lagrange_basis(edgeElemType,pt);  % element shape functions
        J0=dNdxi'*node2(sctr,:);
        detJ0=norm(J0);
        yPt=N'*node2(sctr,2);
        fyPt=-P*(c^2-yPt^2)/(2*I0);
        %f(sctry)=f(sctry)+N*fyPt*detJ0*wt;
        f(sctrx)=f(sctrx)+N*1000*detJ0*wt;
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

scaleFact=100.;
fn=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT DEFORMED DISPLACEMENT PLOT
%Ux = U(xs);
%Uy = U(ys);

figure(fn)
clf
%plot_field(node1+scaleFact*[Ux Uy],element1,elemType,U(ys));
%plot_field(node2+scaleFact*[Ux Uy],element2,elemType,U(ys));

hold on
plot_mesh(node1+scaleFact*[Ux1 Uy1],element1,elemType,'g.-',2.4);
plot_mesh(node2+scaleFact*[Ux2 Uy2],element2,elemType,'r.-',1);
plot_field(node1+scaleFact*[Ux1 Uy1],element1,elemType,Ux1);
plot_field(node2+scaleFact*[Ux2 Uy2],element2,elemType,Ux2);
colorbar
fn=fn+1;
title('DEFORMED DISPLACEMENT IN Y-DIRECTION')

% Comapre numerical displacement to exact value

y = 0;
x = L;
D = 2*c;
uyExact = -P/(6*E0*I0)*(3*nu0*y^2*(L-x)+(4+5*nu0)*D^2*x/4+...
          (3*L-x)*x^2);
      
Uy2(end)
uyExact      

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
plot_mesh(node2+scaleFact*[Ux2 Uy2],element2,elemType,'r.-',1);
hold on
%plot_mesh(node+scaleFact*[U(xs) U(ys)],element,elemType,'g.-');
%plot_mesh(node,element,elemType,'w--');
colorbar
fn=fn+1;
title('DEFORMED STRESS PLOT, BENDING COMPONENT')
