
addpath('~/code/xfem-efg-matlab/fem_util');
addpath ../gmshFiles/
addpath ../post-processing/


clear
colordef black
state = 0;

% ******************************************************************************
% ***                            I N P  U T                                  ***
% ******************************************************************************
tic;
disp('************************************************')
disp('***          S T A R T I N G    R  U N        ***')
disp('************************************************')
disp([num2str(toc),'  START'])

E0          = 1e5;  % Young modulus
nu0         = 0.3;  % Poisson’s ratio
stressState = 'PLANE_STRESS';
a           = 0.5; % hole radius 

plotMesh    = 1;     
computeStr  = 1;
                            
% ******************************************************************************
% ***                    P R E - P R O  C E S S I N G                        ***
% ******************************************************************************


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE ELASTICITY MATRIX
if ( strcmp(stressState,'PLANE_STRESS')  )      % Plane Strain case
  C=E0/(1-nu0^2)*[  1      nu0          0;
                  nu0        1          0;
                    0        0  (1-nu0)/2  ];
else                                            % Plane Strain case
  C=E0/(1+nu0)/(1-2*nu0)*[ 1-nu0      nu0        0;
                            nu0    1-nu0        0;
                              0        0  1/2-nu0 ];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE FINITE ELEMENT MESH
%
meshFile = '11-voided.msh';
mesh     = load_gmsh (meshFile);

elemType = 'Q4';
numnode  = mesh.nbNod;
numelem  = mesh.nbQuads;
node     = mesh.POS(:,1:2);
element  = mesh.QUADS(1:numelem,1:4);

% Finding node groups for boundary conditions

ngr1 = find(mesh.LINES(:,3)==1);
ngr2 = find(mesh.LINES(:,3)==2);
ngr3 = find(mesh.LINES(:,3)==3);
ngr4 = find(mesh.LINES(:,3)==4);

leftNodes  = unique(mesh.LINES(ngr1,1:2)); % nodes
rightNodes = unique(mesh.LINES(ngr3,1:2)); % nodes


vNodes     = [leftNodes(3); rightNodes(3)];
uFixed     = [zeros(length(leftNodes),1);...
              ones(length(rightNodes),1)];
vFixed     = zeros(length(vNodes),1);

%PLOT MESH

if ( plotMesh )  % if plotMesh==1 we will plot the mesh
  clf
  plot_mesh(node,element,elemType,'g.-');
  hold on
  axis off
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DEFINE SYSTEM DATA STRUCTURES

disp([num2str(toc),'  INITIALIZING DATA STRUCTURES'])

U = zeros(2*numnode,1);          % nodal displacement vector
f = zeros(2*numnode,1);          % external  load vector
K = sparse(2*numnode,2*numnode); % stiffness  matrix

% a vector of indicies that quickly address  the x and y portions of the data
% strtuctures so U(xs) returns U_x the nodal  x-displacements

xs=1:numnode;                  % x portion  of u and v vectors
ys=(numnode+1):2*numnode;      % y portion  of u and v vectors

% ******************************************************************************
% ***                          P R O C E  S S I N G                          ***
% ******************************************************************************

%%%%%%%%%%%%%%%%%%%%% COMPUTE STIFFNESS  MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp([num2str(toc),'  COMPUTING STIFFNESS  MATRIX'])

[W,Q]=quadrature(2, 'GAUSS', 2 ); % 2x2 Gaussian quadrature

for e=1:numelem                          % start of element loop
  sctr=element(e,:);          %  element scatter vector
  sctrB=[ sctr sctr+numnode ]; %  vector that scatters a B matrix
  nn=length(sctr);
  B=zeros(3,2*nn);
  for q=1:size(W,1)                        % quadrature loop
    pt=Q(q,:);                              % quadrature point
    wt=W(q);                                % quadrature weight
    [N,dNdxi]=lagrange_basis(elemType,pt);  % element shape functions
    J0=node(sctr,:)'*dNdxi;                % element Jacobian matrix
    invJ0=inv(J0);
    dNdx=dNdxi*invJ0;
    j=det(J0); if j <0 j=-j; end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % COMPUTE B MATRIX
    %        _                                      _
    %        |  N_1,x  N_2,x  ...      0      0  ... |
    %  B  =  |      0      0  ... N_1,y  N_2,y  ... |
    %        |  N_1,y  N_2,y  ... N_1,x  N_2,x  ... |
    %        -                                      -
 
    B(1,1:nn)       = dNdx(:,1)';
    B(2,nn+1:2*nn)  = dNdx(:,2)';
    B(3,1:nn)       = dNdx(:,2)';
    B(3,nn+1:2*nn)  = dNdx(:,1)';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % COMPUTE ELEMENT STIFFNESS AT QUADRATURE  POINT
    K(sctrB,sctrB)=K(sctrB,sctrB)+B'*C*B*W(q)*j;
  end  % of quadrature loop

end    % of element loop
%%%%%%%%%%%%%%%%%%% END OF STIFFNESS  MATRIX COMPUTATION %%%%%%%%%%%%%%%%%%%%%%

% APPLY ESSENTIAL BOUNDARY CONDITIONS
disp([num2str(toc),'  APPLYING BOUNDARY  CONDITIONS'])
bcwt=mean(diag(K)); % a measure of the average  size of an element in K
                    % used to keep the  conditioning of the K matrix
udofs=[leftNodes;rightNodes];          % global indecies  of the fixed x displacements
vdofs=vNodes+numnode;  % global indecies  of the fixed y displacements
f=f-K(:,udofs)*uFixed;  % modify the  force vector
f=f-K(:,vdofs)*vFixed;
f(udofs)=bcwt*uFixed;
f(vdofs)=bcwt*vFixed;
K(udofs,:)=0;  % zero out the rows and  columns of the K matrix
K(vdofs,:)=0;
K(:,udofs)=0;
K(:,vdofs)=0;
K(udofs,udofs)=bcwt*speye(length(udofs));  % put ones*bcwt on the diagonal
K(vdofs,vdofs)=bcwt*speye(length(vdofs));

% SOLVE SYSTEM
disp([num2str(toc),'  SOLVING SYSTEM'])
U=K\f;

%******************************************************************************
%***                    P O S T  -  P R O C E S S I N G                    ***
%******************************************************************************
%
% Here we plot the stresses and displacements  of the solution. As with the
% mesh generation section we don't go  into too much detail - use help
% 'function name' to get more details.
disp([num2str(toc),'  POST-PROCESSING'])

Ux = U(xs);
Uy = U(ys);
scaleFact=0.3;
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

% export the figure to EPS file

% opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',13);
% exportfig(gcf,'plate-q4.eps',opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE STRESS

if (computeStr)
    
    stress=zeros(numelem,size(element,2),3);
    
    stressPoints=[-1 -1;1 -1;1 1;-1  1];
    
    for e=1:numelem                          % start of element loop
        sctr=element(e,:);
        sctrB=[sctr sctr+numnode];
        nn=length(sctr);
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
    stressComp=1;
    figure(fn)
    clf
    plot_field(node+scaleFact*[U(xs) U(ys)],element,elemType,stress(:,:,stressComp));
    hold on
    colorbar
    fn=fn+1;
    title('DEFORMED STRESS PLOT, SIGMA XX')
    
    %exportfig(gcf,'plate-Q4-stress.eps',opts)
end

sigmaXX = zeros(size(node,1),1);
sigmaYY = zeros(size(node,1),1);
sigmaXY = zeros(size(node,1),1);

for e=1:size(element,1)
    connect = element(e,:);
    for in=1:4
        nid = connect(in);
        sigmaXX(nid) = stress(e,in,1);
        sigmaYY(nid) = stress(e,in,2);
        sigmaXY(nid) = stress(e,in,3);
    end
end

VTKPostProcess(node,element,2,'Quad4','../results/plateHoleFullFEM',...
             [sigmaXX sigmaYY sigmaXY],[Ux Uy]);
         
disp([num2str(toc),'  RUN FINISHED'])
% ***************************************************************************
% ***                    E N D  O F    P R O G R A M                    ***
% ***************************************************************************
disp('************************************************')
disp('***            E N D    O F    R U N        ***')
disp('************************************************') 