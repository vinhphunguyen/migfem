% This file implements the Nitsche method to join two 3D meshes.
%
% Discretisation: standard Lagrange finite elements.
%
% Problem: Timoshenko beam in bending.
% Solved by continuous Galerkin to compare with Nitsche.
% 
% Vinh Phu Nguyen
% Cardiff University, UK
% 7 June 2013

addpath ../fem_util/
addpath ../gmshFiles/
addpath ../post-processing/
addpath ../fem-functions/
addpath ../analytical-solutions/
addpath ../meshing/
addpath ../nurbs-util/
addpath ../delamination/

clear all
colordef white
state = 0;
tic;

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,'splinecurve.eps',opts)


% MATERIAL PROPERTIES
E0  = 1000;  % Young modulus
nu0 = 0.3;   % Poisson ratio

% BEAM PROPERTIES
L     = 5;     % length of the beam
width = 1;
t     = 1;     % thicknes

% TIP LOAD
P = -1; 
ubar=-1;

% Constitutive matrices: 
C=zeros(6,6);
C(1:3,1:3)=E0/(1+nu0)/(1-2*nu0)*[ 1-nu0 nu0 nu0;
                                  nu0 1-nu0 nu0;
                                  nu0 nu0 1-nu0];
C(4:6,4:6)=E0/2/(1+nu0)*eye(3);

noGPs = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GENERATE FINITE ELEMENT MESH
%
plotMesh  = 1;
disp([num2str(toc),'   GENERATING MESH'])

% MESH PROPERTIES

elemType  = 'B8'; % the element type for solid1;
elemTypeB = 'Q4'; % the element type for surface;

% domain 1 ----------------------------------------------------------------

controlPts          = zeros(4,2,2,2);

controlPts(1:3,1,1,1) = [0;0;0];
controlPts(1:3,2,1,1) = [L;0;0];
controlPts(1:3,1,2,1) = [0;width;0];
controlPts(1:3,2,2,1) = [L;width;0];

controlPts(1:3,1,1,2) = [0;0;t];
controlPts(1:3,2,1,2) = [L;0;t];
controlPts(1:3,1,2,2) = [0;width;t];
controlPts(1:3,2,2,2) = [L;width;t];

controlPts(4,:,:)   = 1;

uKnot = [0 0 1 1];
vKnot = [0 0 1 1];
wKnot = [0 0 1 1];

solid = nrbmak(controlPts,{uKnot vKnot wKnot});

% evaluate order

%solid = nrbdegelev(solid,[0 0]);


uKnot     = cell2mat(solid.knots(1));
vKnot     = cell2mat(solid.knots(2));
wKnot     = cell2mat(solid.knots(3));

refineCountX = 6;
for i=1:refineCountX
    uKnotVectorU = unique(uKnot);
    uKnotVectorV = unique(vKnot);
    uKnotVectorW = unique(wKnot);
    newKnotsY = [];
    newKnotsX = uKnotVectorU(1:end-1) + 0.5*diff(uKnotVectorU);
    newKnotsY = uKnotVectorV(1:end-1) + 0.5*diff(uKnotVectorV);
    newKnots  = {newKnotsX [] []};
    solid     = nrbkntins(solid,newKnots);
    uKnot     = cell2mat(solid.knots(1));
    vKnot     = cell2mat(solid.knots(2));
    wKnot     = cell2mat(solid.knots(3));
end
newKnots  = {[] [0.5] [0.25 0.5 0.75]};
solid     = nrbkntins(solid,newKnots);

mesh      = buildIGA3DMesh (solid);

node1     = mesh.controlPts;
element1  = mesh.locElems;
element1(:,[3 4 7 8]) = element1(:,[4 3 8 7]);

numx1 = mesh.noElemsU;
numy1 = mesh.noElemsV;
numz1 = mesh.noElemsW;


%% ------------------------------------------------------------------------
% boundary nodes

fixedNode = find(node1(:,1)==0);  
rightNode = find(node1(:,1)==L); 

uNode = [fixedNode; rightNode];

uFixed    = zeros(1,length(fixedNode));  % a vector of u_x for the nodes
vFixed    = zeros(1,length(fixedNode));     
wFixed    = [zeros(1,length(fixedNode)) ubar*ones(1,length(fixedNode))]; 

udofs     = 3*fixedNode-2;
vdofs     = 3*fixedNode-1;
wdofs     = 3*uNode-0;

numnodes = size(node1,1);

% PLOT MESH
if ( plotMesh )  
    clf
    hold on
    plot_mesh(node1,element1,elemType,'g.-',2.4);    
    plot3(node1(fixedNode,1),node1(fixedNode,2),node1(fixedNode,3),'r*');
%    plot3(node2(rightNode,1),node2(rightNode,2),node2(rightNode,3),'rs');    
    axis on
    %axis([0 L -c c])
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE SYSTEM DATA STRUCTURES

f=zeros(3*numnodes,1);          % external load vector
K=zeros(3*numnodes,3*numnodes);  % stiffness matrix

% ******************************************************************************
% ***                          P R O C E S S I N G                           ***
% ******************************************************************************
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%% COMPUTE STIFFNESS MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp([num2str(toc),'   COMPUTING STIFFNESS MATRIX'])

[W,Q]=quadrature( 2, 'GAUSS', 3 ); % 2x2x2 Gaussian quadrature

%% for domain 1 (continuum)

for e=1:size(element1,1)                          % start of element loop
    sctr            = element1(e,:);              % element scatter vector
    nn              = length(sctr);
    sctrB(1:3:3*nn) = 3*sctr-2;
    sctrB(2:3:3*nn) = 3*sctr-1;
    sctrB(3:3:3*nn) = 3*sctr-0;
    pts             = node1(sctr,:);
    for q=1:size(W,1)
        pt=Q(q,:);
        wt=W(q);
        [N,dNdxi]=lagrange_basis(elemType,pt); % element shape functions
        J0    = pts'*dNdxi;                % element Jacobian matrix
        invJ0 = inv(J0);
        dNdx  = dNdxi*invJ0;        
        B     = getBmatrix3D(nn,dNdx);                
        K(sctrB,sctrB)=K(sctrB,sctrB)+B'*C*B*W(q)*det(J0);
    end  % of quadrature loop
end    

%%%%%%%%%%%%%%%%%%% END OF STIFFNESS MATRIX COMPUTATION %%%%%%%%%%%%%%%%%%%


% APPLY ESSENTIAL BOUNDARY CONDITIONS
disp([num2str(toc),'   APPLYING BOUNDARY CONDITIONS'])

[K,f]=applyDirichletBCs3D(K,f,udofs,vdofs,wdofs,uFixed,vFixed,wFixed);

% SOLVE SYSTEM
disp([num2str(toc),'   SOLVING SYSTEM'])
U=K\f;

%******************************************************************************
%*** POST - PROCESSING *** 
%***************************************************

numnode1 = size(node1,1);

index1   = 1:numnode1;

Ux1      = U(3*index1-2);
Uy1      = U(3*index1-1);
Uz1      = U(3*index1-0);


% Here we plot the stresses and displacements of the solution.
disp([num2str(toc),'   POST-PROCESSING'])

dispNorm=L/max(sqrt(Ux1.^2+Uy1.^2));
scaleFact=0.05*dispNorm;

colordef white
figure
clf
hold on
%plot_field(node1+scaleFact*[Ux1 Uy1 Uz1],element1,elemType,Uy1);
%plot_field(node2+scaleFact*[Ux2 Uy2 Uz2],element2,elemType,Uy2);
plot_mesh(node1+scaleFact*[Ux1 Uy1 Uz1],element1,elemType,'r-',1);
%colorbar
axis off
%title('DEFORMED DISPLACEMENT IN Y-DIRECTION')

%% 

stress=zeros(size(element1,1),size(element1,2),6);

stressPoints=[-1 -1 -1;1 -1 -1;1 1 -1;-1 1 -1; ...
              -1 -1  1;1 -1  1;1 1 1; -1 1 1];

for e=1:size(element1,1)
    sctr=element1(e,:);
    nn=length(sctr);
    sctrB(1:3:3*nn) = 3*sctr-2;
    sctrB(2:3:3*nn) = 3*sctr-1;
    sctrB(3:3:3*nn) = 3*sctr-0;    
    for q=1:nn
        pt=stressPoints(q,:);
        [N,dNdxi]=lagrange_basis(elemType,pt);
        J0=node1(sctr,:)'*dNdxi;
        invJ0=inv(J0);
        dNdx=dNdxi*invJ0;
        B      = getBmatrix3D(nn,dNdx);  
        strain=B*U(sctrB);
        stress(e,q,:)=C*strain;
    end
end   % of element loop

vtuFile1    = '../results/fem-3D-3D';

numNode  = size(node1,1);

% normal stresses
sigmaXX = zeros(numNode,2);
sigmaYY = zeros(numNode,2);
sigmaZZ = zeros(numNode,2);

% shear stresses
sigmaXY = zeros(numNode,2);
sigmaYZ = zeros(numNode,2);
sigmaZX = zeros(numNode,2);
% von Mises stress
sigmaVM = zeros(numNode,1);

for e=1:size(element1,1)
    connect = element1(e,:);
    for in=1:8
        nid = connect(in);
        sigmaXX(nid,:) = sigmaXX(nid,:) + [stress(e,in,1) 1];
        sigmaYY(nid,:) = sigmaYY(nid,:) + [stress(e,in,2) 1];
        sigmaZZ(nid,:) = sigmaZZ(nid,:) + [stress(e,in,3) 1];
        sigmaXY(nid,:) = sigmaXY(nid,:) + [stress(e,in,4) 1];
        sigmaYZ(nid,:) = sigmaYZ(nid,:) + [stress(e,in,5) 1];
        sigmaZX(nid,:) = sigmaZX(nid,:) + [stress(e,in,6) 1];
    end
end

% Average nodal stress values (learned from Mathiew Pais XFEM code)
sigmaXX(:,1) = sigmaXX(:,1)./sigmaXX(:,2); sigmaXX(:,2) = [];
sigmaYY(:,1) = sigmaYY(:,1)./sigmaYY(:,2); sigmaYY(:,2) = [];
sigmaZZ(:,1) = sigmaZZ(:,1)./sigmaZZ(:,2); sigmaZZ(:,2) = [];
sigmaXY(:,1) = sigmaXY(:,1)./sigmaXY(:,2); sigmaXY(:,2) = [];
sigmaYZ(:,1) = sigmaYZ(:,1)./sigmaYZ(:,2); sigmaYZ(:,2) = [];
sigmaZX(:,1) = sigmaZX(:,1)./sigmaZX(:,2); sigmaZX(:,2) = [];

VTKPostProcess3d(node1,element1,'B8',vtuFile1,...
    [sigmaXX sigmaYY sigmaZZ sigmaXY sigmaYZ sigmaZX sigmaVM],[Ux1 Uy1 Uz1]);




