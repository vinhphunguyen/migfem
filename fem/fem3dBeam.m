%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Isogeometric analysis for three dimensional elasticity problems.
%
% 3D beam in bending solved with FEM
%
% Vinh Phu Nguyen,
% Johns Hopkins University
% nvinhphu@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('~/code/xfem-efg-matlab/fem_util');
addpath ../gmshFiles/
addpath ../post-processing/

clc
clear all

E0           = 1e5;  % Young modulus
nu0          = 0.3;  % Poisson’s ratio


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE COMPLIANCE MATRIX
C=zeros(6,6);
C(1:3,1:3)=E0/(1+nu0)/(1-2*nu0)*[ 1-nu0 nu0 nu0;
                                  nu0 1-nu0 nu0;
                                  nu0 nu0 1-nu0];
C(4:6,4:6)=E0/(1+nu0)*eye(3);

tic;

a=10;
b=4;
c=2;

nnx = 21;
nny = 7;
nnz = 3;

[node,element]=makeB8mesh(a,b,c,nnx,nny,nnz);

noNodes = size(node,1);
noDofs  = size(node,1) * 3;
noElems = size(element,1);

% Boundary nodes
% find boundary nodes for bounjdary conditions

leftNodes   = find(node(:,1)==0)';
rightNodes  = find(node(:,1)==a)';

% essential boundary conditions

uFixed     = zeros(size(leftNodes));
vFixed     = zeros(size(leftNodes));
wFixed     = zeros(size(leftNodes));


% initialization

K = sparse(noDofs,noDofs);  % global stiffness matrix
u = zeros(noDofs,1);        % displacement vector
f = zeros(noDofs,1);        % external force vector

%jacob = zeros(2,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gauss quadrature rule
[W,Q]=quadrature(  3, 'GAUSS', 3 ); 

% Assembling system of equation
% Stiffness matrix and external force vector

disp([num2str(toc),'  ASSEMBLING THE SYSTEM'])

% Loop over elements (knot spans)

for e=1:noElems
    sctr   = element(e,:);  %  element scatter vector
    sctrB  = [sctr sctr+noNodes sctr+2*noNodes]; % scatters a B matrix
    nn     = length(sctr);
    
    B      = zeros(6,3*nn);
    
    % loop over Gauss points
    
    for gp=1:size(W,1)
        pt      = Q(gp,:);
        wt      = W(gp);
    
        [N, dNdxi] = lagrange_basis('B8',pt); 
        J0         = node(sctr,:)'*dNdxi;
        invJ0      = inv(J0);
        dNdx       = dNdxi*invJ0;
        detJ0      = det(J0);
        
        % B matrix
        
        B(1,1:nn)         = dNdx(:,1)';
        B(2,nn+1:2*nn)    = dNdx(:,2)';
        B(3,2*nn+1:3*nn)  = dNdx(:,3)';
        
        B(4,1:nn)         = dNdx(:,2)';
        B(4,nn+1:2*nn)    = dNdx(:,1)';
        
        B(5,2*nn+1:3*nn)  = dNdx(:,2)';
        B(5,nn+1:2*nn)    = dNdx(:,3)';
        
        B(6,1:nn)         = dNdx(:,3)';
        B(6,2*nn+1:3*nn)  = dNdx(:,1)';
        
        % compute elementary stiffness matrix and
        % assemble it to the global matrix
        
        K(sctrB,sctrB) = K(sctrB,sctrB) + B' * C * B * detJ0 * wt;
    end
end

% Computing external force

f(rightNodes+2*noNodes)=-100;

disp([num2str(toc),'  APPLYING BOUNDARY  CONDITIONS'])

udofs = leftNodes;             % global indecies  of the fixed x disps
vdofs = leftNodes+noNodes;    % global indecies  of the fixed y disps
wdofs = leftNodes+2*noNodes;  % global indecies  of the fixed z disps

bcwt=mean(diag(K)); % a measure of the average  size of an element in K
% used to keep the  conditioning of the K matrix

f=f-K(:,udofs)*uFixed';  % modify the  force vector
f=f-K(:,vdofs)*vFixed';
f=f-K(:,wdofs)*wFixed';

f(udofs) = bcwt*uFixed;
f(vdofs) = bcwt*vFixed;
f(wdofs) = bcwt*wFixed;

K(udofs,:)=0;  % zero out the rows and  columns of the K matrix
K(vdofs,:)=0;
K(wdofs,:)=0;

K(:,udofs)=0;
K(:,vdofs)=0;
K(:,wdofs)=0;

K(udofs,udofs)=bcwt*speye(length(udofs));  % put ones*bcwt on the diagonal
K(vdofs,vdofs)=bcwt*speye(length(vdofs));
K(wdofs,wdofs)=bcwt*speye(length(wdofs));

% SOLVE SYSTEM

disp([num2str(toc),'  SOLVING THE SYSTEM'])
U=K\f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% POST-PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ux    = U(1:noNodes);
Uy    = U(noNodes+1:2*noNodes);
Uz    = U(2*noNodes+1:noDofs);


plot_mesh(node,element,'B8','g.-');
%plot_field(node,elementV,'B8',Uz); do not support B8 elements!!!
view(3)

sigmaXX = zeros(size(node,1),1);
sigmaYY = zeros(size(node,1),1);
sigmaXY = zeros(size(node,1),1);

VTKPostProcess3d(node,element,'B8','../results/beam3DFEM',...
             [sigmaXX sigmaYY sigmaXY sigmaXX sigmaXX sigmaXX],[Ux Uy Uz]);



