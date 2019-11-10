addpath ../fem_util/
addpath ../gmshFiles/
addpath ../post-processing/
addpath ../fem-functions/
addpath ../analytical-solutions/

clear all
colordef black
state = 0;

% ******************************************************************************
% ***                             I N P U T                                  ***
% ******************************************************************************
tic;


% MATERIAL PROPERTIES
E0  = 3e7;  % Young?s modulus
nu0 = 0.3;  % Poisson?s ratio
% BEAM PROPERTIES
L  = 48;     % length of the beam
c  = 3;      % the distance of the outer fiber of the beam from the mid-line

% MESH PROPERTIES
elemType = 'Q4'; % the element type used in the FEM simulation;
numy     = 10;
numx     = 160;
plotMesh = 1;

% TIP LOAD
P = 1000; % the peak magnitude of the traction at the right edge
% STRESS ASSUMPTION


% ******************************************************************************
% ***                    P R E - P R O C E S S I N G                         ***
% ******************************************************************************
I0=2*c^3/3;     % the second polar moment of inertia of the beam cross-section.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stressState = 'PLANE_STRESS';
% COMPUTE ELASTICITY MATRIX
C = elasticityMatrix(E0,nu0,stressState);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE FINITE ELEMENT MESH
% These connectivity matricies refer to the node numbers defined in the
% coordinate matrix node.
disp([num2str(toc),'   GENERATING MESH'])
switch elemType
    case 'Q4'           % here we generate the mesh of Q4 elements
        nnx=numx+1;
        nny=numy+1;
        node=square_node_array([0 -c],[L -c],[L c],[0 c],nnx,nny);
        inc_u=1;
        inc_v=nnx;
        node_pattern=[ 1 2 nnx+2 nnx+1 ];
        element=make_elem(node_pattern,numx,numy,inc_u,inc_v);
    case 'Q9'           % here we generate a mehs of Q9 elements
        nnx=2*numx+1;
        nny=2*numy+1;
        node=square_node_array([0 -c],[L -c],[L c],[0 c],nnx,nny);
        inc_u=2;
        inc_v=2*nnx;
        node_pattern=[ 1 3 2*nnx+3 2*nnx+1 2 nnx+3 2*nnx+2 nnx+1 nnx+2 ];
        element=make_elem(node_pattern,numx,numy,inc_u,inc_v);
    otherwise %?T3?    % and last but not least T3 elements
        nnx=numx+1;
        nny=numy+1;
        node=square_node_array([0 -c],[L -c],[L c],[0 c],nnx,nny);
        node_pattern1=[ 1 2 nnx+1 ];
        node_pattern2=[ 2 nnx+2 nnx+1 ];
        inc_u=1;
        inc_v=nnx;
        element=[make_elem(node_pattern1,numx,numy,inc_u,inc_v);
            make_elem(node_pattern2,numx,numy,inc_u,inc_v) ];
end
% DEFINE BOUNDARIES

%   Here we define the boundary discretizations.
uln=nnx*(nny-1)+1;
urn=nnx*nny;
lrn=nnx;
lln=1;
cln=nnx*(nny-1)/2+1;  % node number at (0,0)
switch elemType
    
    % upper left node number
    % upper right node number
    % lower right node number
    % lower left node number
    case 'Q9'
        rightEdge=[ lrn:2*nnx:(uln-1); (lrn+2*nnx):2*nnx:urn; (lrn+nnx):2*nnx:urn ]';
        leftEdge =[ uln:-2*nnx:(lrn+1); (uln-2*nnx):-2*nnx:1; (uln-nnx):-2*nnx:1 ]';
        edgeElemType='L3';
    otherwise  % same discretizations for Q4 and T3 meshes
        rightEdge=[ lrn:nnx:(uln-1); (lrn+nnx):nnx:urn ]';
        leftEdge =[ uln:-nnx:(lrn+1); (uln-nnx):-nnx:1 ]';
        edgeElemType='L2';
end
% GET NODES ON DISPLACEMENT BOUNDARY
%      Here we get the nodes on the essential boundaries

fixedNode=unique(leftEdge);  % a vector of the node numbers which are fixed in
% the x direction
rightNode=unique(rightEdge);  % a vector of the node numbers which are fixed in

uFixed=zeros(1,length(fixedNode))';     % a vector of the x-displacement for the nodes
vFixed=zeros(1,length(fixedNode))';     % and the y-displacements for fixedNodeY

% for i=1:length(fixedNode)
%     inode = fixedNode(i);
%     pts   = node(inode,:);
%     ux    = P*pts(2)/(6*E0*I0)*(2+nu0)*(pts(2)^2-c^2);
%     uy    = -P/(2*E0*I0)*nu0*pts(2)^2*L;
%     uFixed(i)=ux;
%     vFixed(i)=uy;
% end

numnode=size(node,1);    % number of nodes
numelem=size(element,1); % number of elements
% PLOT MESH
if ( plotMesh )  % if plotMesh==1 we will plot the mesh
    clf
    plot_mesh(node,element,elemType,'g.-',1);
    hold on
    plot_mesh(node,rightEdge,edgeElemType,'bo-',1);
    plot_mesh(node,leftEdge,edgeElemType,'bo-',1);
    plot(node(fixedNode,1),node(fixedNode,2),'r>');
    axis off
    axis([0 L -c c])
end

midNodes  = find(abs(node(:,2))<1e-14);  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE SYSTEM DATA STRUCTURES
%
% Here we define the system data structures
%   f - is the nodal force vector.  It?s structure is the same as U,
%
% K-
%
disp([num2str(toc),' INITIALIZING DATA STRUCTURES'])
U=zeros(2*numnode,1); % nodal displacement vector
f=zeros(2*numnode,1);          % external load vector
K=sparse(2*numnode,2*numnode); % stiffness matrix
% a vector of indicies that quickly address the x and y portions of the data
% strtuctures so U(xs) returns U_x the nodal x-displacements
xs=1:numnode;                  % x portion of u and v vectors
ys=(numnode+1):2*numnode;      % y portion of u and v vectors
% ******************************************************************************
% ***                          P R O C E S S I N G                           ***
% ******************************************************************************
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE EXTERNAL FORCES
%   integrate the tractions on the left and right edges
disp([num2str(toc),'   COMPUTING EXTERNAL LOADS'])
switch elemType  % define quadrature rule
    case 'Q9'
        [W,Q]=quadrature( 4, 'GAUSS', 1 ); %  four point quadrature
    otherwise
        [W,Q]=quadrature( 3, 'GAUSS', 1 ); % three point quadrature
end

% RIGHT EDGE
% for e=1:size(rightEdge,1) % loop over the elements in the right edge
%     sctr=rightEdge(e,:);  % scatter vector for the element
%     sctrx=sctr;           % x scatter vector
%     sctry=sctrx+numnode;  % y scatter vector
%     for q=1:size(W,1)
%         pt=Q(q,:);
%         wt=W(q);
%         [N,dNdxi]=lagrange_basis(edgeElemType,pt);  % element shape functions
%         J0=dNdxi'*node(sctr,:);
%         detJ0=norm(J0);
%         % element Jacobian
%         % determiniat of jacobian
%         yPt=N'*node(sctr,2);
%         fyPt=-P*(c^2-yPt^2)/(2*I0);
%         fyPt=-P;
%         f(sctry)=f(sctry)+N*fyPt*detJ0*wt;  %tor        
%     end % of quadrature loop
% end  % of element loop

f(2*rightNode(end))=-P;

%%%%%%%%%%%%%%%%%%%%% COMPUTE STIFFNESS MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp([num2str(toc),'   COMPUTING STIFFNESS MATRIX'])
switch elemType  % define quadrature rule
    case 'Q9'
        [W,Q]=quadrature( 4, 'GAUSS', 2 ); % 4x4 Gaussian quadrature
    case 'Q4'
        [W,Q]=quadrature( 2, 'GAUSS', 2 ); % 2x2 Gaussian quadrature
    otherwise
        [W,Q]=quadrature( 1, 'TRIANGULAR', 2 ); % 1 point triangural quadrature
end
for e=1:numelem                          % start of element loop
    sctr=element(e,:);           % element scatter vector
    sctrB=[ sctr sctr+numnode ]; % vector that scatters a B matrix
    nn=length(sctr);
    for q=1:size(W,1)
        pt=Q(q,:);
        wt=W(q);
        [N,dNdxi]=lagrange_basis(elemType,pt); % element shape functions
        J0=node(sctr,:)'*dNdxi;                % element Jacobian matrix
        invJ0=inv(J0);
        dNdx=dNdxi*invJ0;
        
        B(1,1:nn)       = dNdx(:,1)';
        B(2,nn+1:2*nn)  = dNdx(:,2)';
        B(3,1:nn)       = dNdx(:,2)';
        B(3,nn+1:2*nn)  = dNdx(:,1)';
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % COMPUTE ELEMENT STIFFNESS AT QUADRATURE POINT
        K(sctrB,sctrB)=K(sctrB,sctrB)+B'*C*B*W(q)*det(J0);
    end  % of quadrature loop
end    % of element loop
%%%%%%%%%%%%%%%%%%% END OF STIFFNESS MATRIX COMPUTATION %%%%%%%%%%%%%%%%%%%%%%

% APPLY ESSENTIAL BOUNDARY CONDITIONS
disp([num2str(toc),'   APPLYING BOUNDARY CONDITIONS'])
bcwt=mean(diag(K)); % a measure of the average size of an element in K
% used to keep the conditioning of the K matrix
udofs=fixedNode;           % global indecies of the fixed x displacements
vdofs=fixedNode+numnode;   % global indecies of the fixed y displacements
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
%*** POST - PROCESSING *** %***************************************************

ptId  = intersect (find(node(:,1)==48),find(node(:,2)==0));
numUy = U(ptId+numnode)


data = [node(midNodes,1) U(midNodes+numnode)];
csvwrite('beam-2D.csv',data);

% Here we plot the stresses and displacements of the solution. As with the
% mesh generation section we don?t go into too much detail - use help
% ?function name? to get more details.
disp([num2str(toc),'   POST-PROCESSING'])
dispNorm=L/max(sqrt(U(xs).^2+U(ys).^2));
scaleFact=0.1*dispNorm;
fn=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT DEFORMED DISPLACEMENT PLOT
Ux = U(xs);
Uy = U(ys);

figure(fn)
clf
plot_field(node+scaleFact*[Ux Uy],element,elemType,U(xs));
hold on
plot_mesh(node+scaleFact*[Ux Uy],element,elemType,'g.-',1);
%plot_mesh(node,element,elemType,'w--');
colorbar
fn=fn+1;
title('DEFORMED DISPLACEMENT IN Y-DIRECTION')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE STRESS
stress=zeros(numelem,size(element,2),4);
switch elemType  % define quadrature rule
    case 'Q9'
        stressPoints=[-1 -1;1 -1;1 1;-1 1;0 -1;1 0;0 1;-1 0;0 0 ];
    case 'Q4'
        stressPoints=[-1 -1;1 -1;1 1;-1 1];
    otherwise
        stressPoints=[0 0;1 0;0 1];
end

for e=1:numelem
    sctr=element(e,:);
    sctrB=[sctr sctr+numnode];
    nn=length(sctr);
    for q=1:nn
        pt=stressPoints(q,:);
        [N,dNdxi]=lagrange_basis(elemType,pt);
        J0=node(sctr,:)'*dNdxi;
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
        sigma =C*strain;
        stress(e,q,1:3)=sigma;
        stress(e,q,4)  = sqrt(sigma(1)^2+sigma(2)^2-sigma(1)*sigma(2)+3*(sigma(3)^2));
    end
end   % of element loop


%print(gcf, '-depsc', '-painters', 'output.eps')
%print(fn,?-djpeg90?,[?beam_?,elemType,?_sigma?,num2str(stressComp),?.jpg?])

numNode  = size(node,1);

% normal stresses
sigmaXX = zeros(numNode,2);
sigmaYY = zeros(numNode,2);
sigmaXY = zeros(numNode,2);
sigmaVM = zeros(numNode,2);

for e=1:size(element,1)
    connect = element(e,:);
    for in=1:4
        nid = connect(in);
        sigmaXX(nid,:) = sigmaXX(nid,:) + [stress(e,in,1) 1];
        sigmaYY(nid,:) = sigmaYY(nid,:) + [stress(e,in,2) 1];
        sigmaXY(nid,:) = sigmaXY(nid,:) + [stress(e,in,3) 1];
        sigmaVM(nid,:) = sigmaVM(nid,:) + [stress(e,in,4) 1];
    end
end

% Average nodal stress values (learned from Mathiew Pais XFEM code)
sigmaXX(:,1) = sigmaXX(:,1)./sigmaXX(:,2); sigmaXX(:,2) = [];
sigmaYY(:,1) = sigmaYY(:,1)./sigmaYY(:,2); sigmaYY(:,2) = [];
sigmaXY(:,1) = sigmaXY(:,1)./sigmaXY(:,2); sigmaXY(:,2) = [];
sigmaVM(:,1) = sigmaVM(:,1)./sigmaVM(:,2); sigmaVM(:,2) = [];

stressComp=3;
set(gca,'FontSize',14)
clf
plot_field(node+scaleFact*[U(xs) U(ys)],element,elemType,sigmaVM);
hold on
%plot_mesh(node+scaleFact*[U(xs) U(ys)],element,elemType,'g.-');
colorbar
axis off
title('\sigma_{xx}')

%% 

aa=40;
elems = aa:numx:aa+numx*(numy-1);

disp    = zeros(length(elems),2);
sigma   = zeros(length(elems),2);
sigmaRef = zeros(length(elems),2);
xcoord     = zeros(length(elems),2);

for e=1:length(elems)
    ie = elems(e);
    sctr=element(ie,:);
    sctrB=[sctr sctr+numnode];
    nn=length(sctr);
    
    pt=[0 0];
    [N,dNdxi]=lagrange_basis(elemType,pt);
    J0=node(sctr,:)'*dNdxi;
    invJ0=inv(J0);
    dNdx=dNdxi*invJ0;
    yPt=N'*node(sctr,:);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % COMPUTE B MATRIX
        B=zeros(3,2*nn);
        B(1,1:nn)       = dNdx(:,1)';
        B(2,nn+1:2*nn)  = dNdx(:,2)';
        B(3,1:nn)       = dNdx(:,2)';
        B(3,nn+1:2*nn)  = dNdx(:,1)';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % COMPUTE ELEMENT STRAIN AND STRESS AT STRESS POINT
        disp(e,1) = N'*U(sctr);
    disp(e,2) = N'*U(sctr+numnode);
    strain=B*U(sctrB);
    stress=C*strain;
    sigma(e,1)    = stress(1);
    sigma(e,2)    = stress(3);
%     sigmaRef(e,1) = 1000/I0*(L-yPt(1))*yPt(2);
%     sigmaRef(e,2) = -1000/2/I0*(c^2-yPt(2)^2);
    xcoord(e,:)     = yPt;
end   % of element loop

save('../continuum-beam/fem.mat','sigma','xcoord','disp');

aa=80;
elems = aa:numx:aa+numx*(numy-1);

sigma   = zeros(length(elems),2);
disp    = zeros(length(elems),2);
sigmaRef = zeros(length(elems),2);
xcoord     = zeros(length(elems),2);

for e=1:length(elems)
    ie = elems(e);
    sctr=element(ie,:);
    sctrB=[sctr sctr+numnode];
    nn=length(sctr);
    
    pt=[-1 -1];
    [N,dNdxi]=lagrange_basis(elemType,pt);
    J0=node(sctr,:)'*dNdxi;
    invJ0=inv(J0);
    dNdx=dNdxi*invJ0;
    yPt=N'*node(sctr,:);
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
    stress=C*strain;
    disp(e,1) = N'*U(sctr);
    disp(e,2) = N'*U(sctr+numnode);
    sigma(e,1)    = stress(1);
    sigma(e,2)    = stress(3);
%     sigmaRef(e,1) = 1000/I0*(L-yPt(1))*yPt(2);
%     sigmaRef(e,2) = -1000/2/I0*(c^2-yPt(2)^2);
    xcoord(e,:)     = yPt;
end   % of element loop

save('../continuum-beam/fem1.mat','sigma','xcoord','disp');

% colordef white
% figure,set (gcf,'Color','w')
% set(gca,'FontSize',14)
% hold on
% plot(xcoord(:,1),sigmaRef(:,1),'k-','LineWidth',1.4);
% plot(xcoord(:,1),sigma(:,1),'o','MarkerEdgeColor','k',...
%     'MarkerFaceColor','g',...
%     'MarkerSize',6.5);
% plot(xcoord(:,1),sigmaRef(:,2),'k-','LineWidth',1.4);
% plot(xcoord(:,1),sigma(:,2),'s','MarkerEdgeColor','k',...
%     'MarkerFaceColor','g',...
%     'MarkerSize',6.5);
% h=legend('exact','coupling');
% xlabel('x')
% ylabel('\sigma_{xx}')
% grid on
%axis([0 5 -0.55 0])

