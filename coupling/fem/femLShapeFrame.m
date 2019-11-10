

addpath ../../fem_util/

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

% MATERIAL PROPERTIES

E0  = 1000;  % Young?s modulus
nu0 = 0.3;  % Poisson?s ratio

stressState ='PLANE_STRESS'; % set  to either 'PLANE_STRAIN' or "PLANE_STRESS'
plotMesh    = 1;
computeStr  = 1;

F           = 1/2; % 6.25 [N/mm]

% ******************************************************************************
% ***                    P R E - P R O  C E S S I N G                        ***
% ******************************************************************************


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE ELASTICITY MATRIX

C = elasticityMatrix(E0,nu0,stressState);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE FINITE ELEMENT MESH
%
meshFile = 'lshape.msh';
mesh     = load_gmsh (meshFile);

elemType = 'Q4';
numnode  = mesh.nbNod;
numelem  = mesh.nbQuads;
node     = mesh.POS(:,1:2);
element  = mesh.QUADS(1:numelem,1:4);

% Finding node groups for boundary conditions

ngr1 = find(mesh.LINES(:,3)==333);
ngr2 = find(mesh.LINES(:,3)==444);

fixedNodes = unique(mesh.LINES(ngr1,1:2)); % nodes
forceEdge  = mesh.LINES(ngr2,1:2);

fixedNodesU = [fixedNodes; unique(forceEdge)];
fixedNodesV = fixedNodes;

uFixed     = zeros(length(fixedNodesU),1);
vFixed     = zeros(length(fixedNodesV),1);

%PLOT MESH

if ( plotMesh )  % if plotMesh==1 we will plot the mesh
    clf
    plot_mesh(node,element,elemType,'g.-',1);
    hold on
    axis off
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DEFINE SYSTEM DATA STRUCTURES

f = zeros(2*numnode,1);          % external  load vector
K = sparse(2*numnode,2*numnode); % stiffness  matrix

% ******************************************************************************
% ***                          P R O C E  S S I N G                          ***
% ******************************************************************************

%%%%%%%%%%%%%%%%%%%%% COMPUTE STIFFNESS  MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp([num2str(toc),'  COMPUTING STIFFNESS  MATRIX'])

[W,Q]=quadrature(  2, 'GAUSS', 2 ); % 2x2 Gaussian quadrature

for e=1:numelem                                 % start of element loop
    sctr=element(e,:);                          %  element scatter vector
    nn=length(sctr);
    sctrB(1:2:2*nn) = 2*sctr-1;
    sctrB(2:2:2*nn) = 2*sctr-0;
    for q=1:size(W,1)                           % quadrature loop
        pt=Q(q,:);                              % quadrature point
        wt=W(q);                                % quadrature weight
        [N,dNdxi]=lagrange_basis(elemType,pt);  % element shape functions
        J0=node(sctr,:)'*dNdxi;                 % element Jacobian matrix
        invJ0=inv(J0);
        dNdx=dNdxi*invJ0;
        B(1,1:2:2*nn)  = dNdx(:,1)';
        B(2,2:2:2*nn)  = dNdx(:,2)';
        B(3,1:2:2*nn)  = dNdx(:,2)';
        B(3,2:2:2*nn)  = dNdx(:,1)';
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % COMPUTE ELEMENT STIFFNESS AT QUADRATURE  POINT
        K(sctrB,sctrB)=K(sctrB,sctrB)+B'*C*B*W(q)*det(J0);
    end  % of quadrature loop
end    % of element loop
%%%%%%%%%%%%%%%%%%% END OF STIFFNESS  MATRIX COMPUTATION %%%%%%%%%%%%%%%%%%%%%%

[W,Q]=quadrature(  2, 'GAUSS', 1 ); % 1 point quadrature

% RIGHT EDGE

for e=1:size(forceEdge,1) % loop over the  elements in the right edge
    sctr  = forceEdge(e,:); % scatter  vector for the element
    sctrx = 2*sctr-1;           % x scatter  vector
    sctry = 2*sctr-0;           % x scatter  vector
    for q=1:size(W,1)                              % quadrature loop
        pt       = Q(q,:);                           % quadrature point
        wt       = W(q);                             % quadrature weight
        [N,dNdxi]=lagrange_basis('L2',pt);           % element shape functions
        J0       = dNdxi'*node(sctr,:);              % element Jacobian
        detJ0    = norm(J0);                         % determiniat of jacobian
        f(sctry) = f(sctry) - N*F*detJ0*wt;  % scatter force into global force vector
    end % of quadrature loop
end  % of element loop


% APPLY ESSENTIAL BOUNDARY CONDITIONS
disp([num2str(toc),'  APPLYING BOUNDARY  CONDITIONS'])
bcwt=mean(diag(K)); % a measure of the average  size of an element in K
% used to keep the  conditioning of the K matrix
udofs=2*fixedNodesU-1;          % global indecies  of the fixed x displacements
vdofs=2*fixedNodesV;  % global indecies  of the fixed y displacements
f=f-K(:,udofs)*uFixed;  % modify the  force vector
f=f-K(:,vdofs)*vFixed;
f(udofs)=uFixed;
f(vdofs)=vFixed;
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

Ux = U(1:2:end);
Uy = U(2:2:end);

scaleFact=10;
fn=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT DEFORMED DISPLACEMENT PLOT
figure(fn)
clf
plot_field(node+scaleFact*[Ux Uy],element,elemType,Ux);
hold on
plot_mesh(node,element,elemType,'g-',1);
colorbar
fn=fn+1;
title('DEFORMED DISPLACEMENT IN X-DIRECTION')

% export the figure to EPS file

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',13);
%exportfig(gcf,'plate-q4.eps',opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE STRESS

if (computeStr)
    
    stress=zeros(numelem,size(element,2),3);
    
    stressPoints=[-1 -1;1 -1;1 1;-1  1];
    
    for e=1:size(element,1)
        sctr=element(e,:);
        sctrB(1:2:2*nn) = 2*sctr-1;
        sctrB(2:2:2*nn) = 2*sctr-0;
        nn=length(sctr);
        Ce=C;
        
        for q=1:nn
            pt=stressPoints(q,:);
            [N,dNdxi]=lagrange_basis(elemType,pt);
            J0=node(sctr,:)'*dNdxi;
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
    colordef white
    stressComp=2;
    figure(fn)
    clf
    plot_field(node+scaleFact*[Ux Uy],element,elemType,stress(:,:,stressComp));
    hold on
    %plot_mesh(node+scaleFact*[Ux Uy],element,elemType,'g-',1);
    colorbar
    axis off
    fn=fn+1;
    title('DEFORMED STRESS PLOT, SIGMA XX')
    
    %exportfig(gcf,'plate-Q4-stress.eps',opts)
end


numNode  = size(node,1);


% normal stresses
sigmaXX = zeros(numNode,2);
sigmaYY = zeros(numNode,2);
% shear stresses
sigmaXY = zeros(numNode,2);


for e=1:size(element,1)
    connect = element(e,:);
    for in=1:4
        nid = connect(in);
        sigmaXX(nid,:) = sigmaXX(nid,:) + [stress(e,in,1) 1];
        sigmaYY(nid,:) = sigmaYY(nid,:) + [stress(e,in,2) 1];
        sigmaXY(nid,:) = sigmaXY(nid,:) + [stress(e,in,3) 1];
    end
end

% Average nodal stress values (learned from Mathiew Pais XFEM code)
sigmaXX(:,1) = sigmaXX(:,1)./sigmaXX(:,2); sigmaXX(:,2) = [];
sigmaYY(:,1) = sigmaYY(:,1)./sigmaYY(:,2); sigmaYY(:,2) = [];
sigmaXY(:,1) = sigmaXY(:,1)./sigmaXY(:,2); sigmaXY(:,2) = [];

figure
plot_field(node+scaleFact*[Ux Uy],element,elemType,sigmaXY);
hold on
%plot_mesh(node+scaleFact*[Ux Uy],element,elemType,'g-',1);
colorbar
axis off

disp([num2str(toc),'  RUN FINISHED'])
% ***************************************************************************
% ***                    E N D  O F    P R O G R A M                    ***
% ***************************************************************************
disp('************************************************')
disp('***            E N D    O F    R U N        ***')
disp('************************************************')


%% save the result to mat file

data.mesh = element;
data.node = node;
data.disp = U;

save('femLShapeFrame.mat', '-struct', 'data');


