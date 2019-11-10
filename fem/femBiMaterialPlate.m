% Finite element code for linear elasticity problems.
% Rectangular plate made of two isotropic materials.
% Only a quater plate is modelled.
% Vinh Phu Nguyen, nvinhphu@gmail.com
% TU Delft, The Netherlands

addpath('~/code/xfem-efg-matlab/fem_util');
addpath ../gmshFiles/
addpath ../nurbs-geopdes/inst/
addpath ../nurbs-util/
addpath ../meshing/
addpath ../fem-functions/
addpath ../post-processing/
addpath ../xiga/

clear all
clc
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


E1           = 1e2;  % Young modulus of matrix
nu1          = 0.2;  % Poisson ratio
E2           = 1e1;  % of inclusion
nu2          = 0.3;
stressState = 'PLANE_STRAIN';

% COMPUTE ELASTICITY MATRIX
C1 = elasticityMatrix(E1,nu1,stressState);
C2 = elasticityMatrix(E2,nu2,stressState);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE FINITE ELEMENT MESH

numx = 62;
numy = 62;
L    = 1;

nnx=numx+1;
nny=numy+1;
node=square_node_array([0 0],[L 0],[L L],[0 L],nnx,nny);
inc_u=1;
inc_v=nnx;
node_pattern=[ 1 2 nnx+2 nnx+1 ];
element=make_elem(node_pattern,numx,numy,inc_u,inc_v);
 
elemType = 'Q4';
numnode  = size(node,1);
numelem  = size(element,1);

% Finding node groups for boundary conditions

bottomNodes = find(node(:,2)==0)';
rightNodes  = find(node(:,1)==L)';
leftNodes   = find(node(:,1)==0)';
topNodes    = find(node(:,2)==L)';

% essential boundary conditions

uNodes = bottomNodes(floor(numx/2)+1);
vNodes = [bottomNodes topNodes];

uFixed     = zeros(size(uNodes));
vFixed     = [zeros(size(bottomNodes)) 0.1*ones(size(topNodes))]';

udofs = uNodes; % global indecies  of the fixed x disps
vdofs = vNodes+numnode;   % global indecies  of the fixed y disps


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

%%%%%%%%%%%%%%%%%%%%% COMPUTE STIFFNESS  MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%
disp([num2str(toc),'  COMPUTING STIFFNESS  MATRIX'])

[W,Q]=quadrature(  2, 'GAUSS', 2 ); % 1 GP rule

for e=1:numelem                          % start of element loop
    sctr=element(e,:);          %  element scatter vector
    sctrB=[ sctr sctr+numnode ]; %  vector that scatters a B matrix
    nn=length(sctr);
    pts=node(sctr,:); 
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
        x         = N' * pts;
        levelset  = x(2) - 0.5;
        
        if levelset >= 1e-10
            C = C1;
        else
            C = C2;
        end
        
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
scaleFact=5;
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

% export the figure to EPS file

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',13);
%exportfig(gcf,'plate-q4.eps',opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE STRESS

if (computeStr)
    
    stress = zeros(numelem,size(element,2),3);
    strain = zeros(numelem,size(element,2),3);
    
    stressPoints=[-1 -1;1 -1;1 1;-1  1]; % Q4 elems
    
    for e=1:numelem                          % start of element loop
        sctr=element(e,:);
        sctrB=[sctr sctr+numnode];
        nn=length(sctr);
        pts=node(sctr,:); 
        
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
            
            x         = N' * pts;
            levelset  = x(2) - 0.5;
            
            if levelset >= 1e-10
                C = C1;
            else
                C = C2;
            end
        
            % COMPUTE ELEMENT STRAIN AND STRESS  AT STRESS POINT
            strainGp     =B*U(sctrB);
            stress(e,q,:)=C*strainGp;
            strain(e,q,:)=strainGp;
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
    
    VTKPostProcess(node,element,2,'Quad4','../results/platInclusionFEM',...
        [sigmaXX sigmaYY sigmaXY],[Ux Uy]);
    
    figure(fn)
    clf
    plot_field(node,element,elemType,sigmaXY);
    plot_mesh(node,element,'Q4','w-');
    hold on
    colorbar
    fn=fn+1;
    title('DEFORMED STRESS IN Y-DIRECTION')
    
    figure(fn)
    clf
    plot_field(node,element,elemType,strain(:,:,3));
    plot_mesh(node,element,'Q4','w-');
    hold on
    colorbar
    fn=fn+1;
    title('DEFORMED STRAIN IN XY-DIRECTION')
end

disp([num2str(toc),'  RUN FINISHED'])
% ***************************************************************************
% ***                    E N D  O F    P R O G R A M                    ***
% ***************************************************************************
disp('************************************************')
disp('***            E N D    O F    R U N        ***')
disp('************************************************')