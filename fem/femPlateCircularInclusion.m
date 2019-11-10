% Finite element code for linear elasticity problems.
% Infinite plate with a centered circular soft inclusion.
% Only a quater plate is modelled.
% Vinh Phu Nguyen, nvinhphu@gmail.com
% TU Delft, The Netherlands

addpath('~/code/xfem-efg-matlab/fem_util');
addpath ../gmshFiles/
addpath ../fem-functions/
addpath ../post-processing/


clear
colordef black
state = 0;

% *************************************************************************
% ***                       I N P  U T                                  ***
% *************************************************************************
tic;
disp('************************************************')
disp('***          S T A R T I N G    R  U N        ***')
disp('************************************************')
disp([num2str(toc),'  START'])

plotMesh    = 1;
computeStr  = 1;
force       = 0; % displacement control

% *************************************************************************
% ***               P R E - P R O  C E S S I N G                        ***
% *************************************************************************


E1  = 1000;  % Young modulus of matrix
E2  = 1;     % Young modulus of inclusion
nu1 = 0.3;  % Poisson ratio
nu2 = 0.3;  % Poisson ratio
stressState ='PLANE_STRAIN'; % either 'PLANE_STRAIN' or "PLANE_STRESS
t0  = -10;

% COMPUTE ELASTICITY MATRIX
C1 = elasticityMatrix(E1,nu1,stressState);
C2 = elasticityMatrix(E2,nu2,stressState);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE FINITE ELEMENT MESH
% Read Gmsh mesh file
meshFile = 'plateInclusion.msh';
mesh     = load_gmsh (meshFile);

elemType = 'T3';
numnode  = mesh.nbNod;
numelem  = mesh.nbTriangles;
node     = mesh.POS(:,1:2);
element  = mesh.TRIANGLES(1:numelem,1:3);
material = mesh.TRIANGLES(1:numelem,4);

% check if Jacobian is negative

element  = tricheck(node,element,1);

% Finding node groups for boundary conditions

ngr1 = find(mesh.LINES(:,3)==111); % fix y disp
ngr2 = find(mesh.LINES(:,3)==222); % fix x disp
ngr3 = find(mesh.LINES(:,3)==333); % force x dir
ngr4 = find(mesh.LINES(:,3)==444); % nothing

fixedXNodes = unique(mesh.LINES(ngr2,1:2)); % nodes
fixedYNodes = unique(mesh.LINES(ngr1,1:2)); % nodes
leftEdge   = mesh.LINES(ngr3,1:2);
topEdge     = mesh.LINES(ngr4,1:2);
leftNodes  = unique(leftEdge);
topNodes    = unique(topEdge);

uFixed     = zeros(length(fixedXNodes),1);
vFixed     = zeros(length(fixedYNodes),1);

udofs=fixedXNodes;          % global indecies  of the fixed x displacements
vdofs=fixedYNodes+numnode;  % global indecies  of the fixed y displacements

if (force==0)
    ubar   = -0.04;
    uFixed = [uFixed; ubar*ones(length(leftNodes),1)];
    udofs  = [udofs; leftNodes];
end

%PLOT MESH

if ( plotMesh )  % if plotMesh==1 we will plot the mesh
    clf
    plot_mesh(node,element,elemType,'g.-');
    hold on
    axis off
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DEFINE SYSTEM DATA STRUCTURES
%
% Here we define the system data structures
%  U - is vector of the nodal displacements it is of length 2*numnode. The
%      displacements in the x-direction are in the top half of U and the
%      y-displacements are in the lower half of U, for example the displacement
%      in the y-direction for node number I is at U(I+numnode)
%  f - is the nodal force vector.  It's structure is the same as U,
%      i.e. f(I+numnode) is the force in the y direction at node I
%  K - is the global stiffness matrix and is structured the same as with U and f
%      so that K_IiJj is at K(I+(i-1)*numnode,J+(j-1)*numnode)

disp([num2str(toc),'  INITIALIZING DATA STRUCTURES'])

U = zeros(2*numnode,1);          % nodal displacement vector
f = zeros(2*numnode,1);          % external  load vector
K = sparse(2*numnode,2*numnode); % stiffness  matrix

% a vector of indicies that quickly address  the x and y portions of the data
% strtuctures so U(xs) returns U_x the nodal  x-displacements

xs=1:numnode;                  % x portion  of u and v vectors
ys=(numnode+1):2*numnode;      % y portion  of u and v vectors

% *************************************************************************
% ***                     P R O C E  S S I N G                          ***
% *************************************************************************

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE EXTERNAL FORCES
%  integrate the external force on the left edge
disp([num2str(toc),'  COMPUTING EXTERNAL  LOADS'])

[W,Q]=quadrature(  2, 'GAUSS', 1 ); % 1 point quadrature

if (force==1)
    for e=1:size(leftEdge,1) % loop over the  elements in the left edge
        sctr  = leftEdge(e,:);  % scatter  vector for the element
        sctrx = sctr;           % x scatter  vector
        sctry = sctr + numnode;
        for q=1:size(W,1)
            pt       = Q(q,:);
            wt       = W(q);
            [N,dNdxi]=lagrange_basis('L2',pt);
            J0       = dNdxi'*node(sctr,:);
            x        = N'*node(sctr,:);
            detJ0    = norm(J0);
            
            f(sctrx) = f(sctrx)+N*t0*detJ0*wt;
        end % of quadrature loop
    end  % of element loop
end


%%%%%%%%%%%%%%%%%%%%% COMPUTE STIFFNESS  MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%
disp([num2str(toc),'  COMPUTING STIFFNESS  MATRIX'])

[W,Q]=quadrature(  1, 'TRIANGULAR', 2 ); % 1 GP rule

for e=1:numelem                          % start of element loop
    sctr=element(e,:);          %  element scatter vector
    sctrB=[ sctr sctr+numnode ]; %  vector that scatters a B matrix
    nn=length(sctr);
    
    % choose correct material
    if material(e)== 10
        C = C1;
    else
        C = C2;
    end
    
    for q=1:size(W,1)                        % quadrature loop
        pt=Q(q,:);                              % quadrature point
        wt=W(q);                                % quadrature weight
        [N,dNdxi]=lagrange_basis(elemType,pt);  % element shape functions
        J0=node(sctr,:)'*dNdxi;                % element Jacobian matrix
        invJ0=inv(J0);
        dNdx=dNdxi*invJ0;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % COMPUTE B MATRIX
        B=zeros(3,2*nn);
        B(1,1:nn)      = dNdx(:,1)';
        B(2,nn+1:2*nn)  = dNdx(:,2)';
        B(3,1:nn)      = dNdx(:,2)';
        B(3,nn+1:2*nn)  = dNdx(:,1)';
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % COMPUTE ELEMENT STIFFNESS AT QUADRATURE  POINT
        K(sctrB,sctrB)=K(sctrB,sctrB)+B'*C*B*W(q)*det(J0);
    end  % of quadrature loop
    
end    % of element loop
%%%%%%%%%%%%%%%%%%% END OF STIFFNESS  MATRIX COMPUTATION %%%%%%%%%%%%%%%%%%%%%%

% APPLY ESSENTIAL BOUNDARY CONDITIONS and SOLVE


solve

%**************************************************************************
%***                 P O S T  -  P R O C E S S I N G                    ***
%**************************************************************************
%
disp([num2str(toc),'  POST-PROCESSING'])

Ux = U(xs);
Uy = U(ys);
scaleFact=10;
fn=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT DEFORMED DISPLACEMENT PLOT
figure(fn)
clf
plot_field(node+scaleFact*[Ux Uy],element,elemType,Ux);
hold on
colorbar
fn=fn+1;
title('DEFORMED DISPLACEMENT IN X-DIRECTION')

figure(fn)
clf
plot_field(node+scaleFact*[Ux Uy],element,elemType,Uy);
hold on
colorbar
fn=fn+1;
title('DEFORMED DISPLACEMENT IN Y-DIRECTION')

uXNode1 = U(4)

% export the figure to EPS file

%opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',13);
%exportfig(gcf,'plate-q4.eps',opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE STRESS

if (computeStr)
    
    stress=zeros(numelem,size(element,2),3);
    
    stressPoints=[-1 -1;1 -1;1 1;-1  1]; % Q4 elems
    stressPoints=[0 0;1 0;0 1]; % T3 elems
    
    for e=1:numelem                          % start of element loop
        sctr=element(e,:);
        sctrB=[sctr sctr+numnode];
        nn=length(sctr);
        % choose correct material
        if material(e)== 10
            C = C1;
        else
            C = C2;
        end
        for q=1:nn
            pt=stressPoints(q,:);                      % stress point
            [N,dNdxi]=lagrange_basis(elemType,pt);    % element shape functions
            J0=node(sctr,:)'*dNdxi;                    % element Jacobian matrix
            invJ0=inv(J0);
            dNdx=dNdxi*invJ0;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % COMPUTE B MATRIX
            B=zeros(3,2*nn);
            B(1,1:nn)      = dNdx(:,1)';
            B(2,nn+1:2*nn)  = dNdx(:,2)';
            B(3,1:nn)      = dNdx(:,2)';
            B(3,nn+1:2*nn)  = dNdx(:,1)';
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % COMPUTE ELEMENT STRAIN AND STRESS  AT STRESS POINT
            strain=B*U(sctrB);
            stress(e,q,:)=C*strain;
        end
    end  % of element loop
    
    % export to VTK format to plot in Mayavi or Paraview
    
    sigmaXX = zeros(size(node,1),1);
    sigmaYY = zeros(size(node,1),1);
    sigmaXY = zeros(size(node,1),1);
    
    for e=1:size(element,1)
        connect = element(e,:);
        for in=1:length(connect)
            nid = connect(in);
            sigmaXX(nid) = stress(e,in,1);
            sigmaYY(nid) = stress(e,in,2);
            sigmaXY(nid) = stress(e,in,3);
        end
    end
    
    VTKPostProcess(node,element,2,'Tri3','../results/platInclusionFEM',...
             [sigmaXX sigmaYY sigmaXY],[Ux Uy]);
end

disp([num2str(toc),'  RUN FINISHED'])
% ***************************************************************************
% ***                    E N D  O F    P R O G R A M                    ***
% ***************************************************************************
disp('************************************************')
disp('***            E N D    O F    R U N        ***')
disp('************************************************')