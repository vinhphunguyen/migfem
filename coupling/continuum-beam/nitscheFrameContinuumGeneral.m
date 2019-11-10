% This file implements the Nitsche method to join two mechanical models:
% a 2D continuum model and a frame model.
% The beam elements can be placed arbitrarily in space.
% This is a prototype for solving frame problems with continuum elements
% in region of high gradients and beam elements in other places.
%
% Dofs: 
%     continuum nodes: ux and uy
%     beam nodes:      ux,uy and rotation
% Assembly: assume that continuum nodes are stored first followed by beam nodes.
%
% Discretisation: standard Lagrange finite elements.
%
% Problem: Timoshenko beam in bending.
% 
% Vinh Phu Nguyen
% Cardiff University, UK
% 7 June 2013

addpath ../../fem_util/
addpath ../../gmshFiles/
addpath ../../post-processing/
addpath ../../fem-functions/
addpath ../../analytical-solutions/


clear all
colordef white
state = 0;
tic;

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,'splinecurve.eps',opts)


% MATERIAL PROPERTIES
E0  = 1000;  % Young?s modulus
nu0 = 0.3;  % Poisson?s ratio

% BEAM PROPERTIES
L  = 8;        % length of the beam
c  = 0.5;      % the distance of the outer fiber of the beam from the mid-line
t  = 2*c;      % thickness, area = bxt
b  = 1; 

I0 = 2*c^3/3;  % the second polar moment of inertia of the beam cross-section.
q  = 0;        % uniformly distributed loads
F  = 1;        % end force
G  = E0/2/(1+nu0);
k  = E0*I0; 
A  = b*t;
shearFactor = 5/6;

plotMesh  = 1;

% TIP LOAD
P = 1000; % the peak magnitude of the traction at the right edge

stressState = 'PLANE_STRESS';
% COMPUTE ELASTICITY MATRIX
C = elasticityMatrix(E0,nu0,stressState);

noGPs = 2;

% penalty parameter

alpha = 1e7;
%alpha=6.1673e+03;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE FINITE ELEMENT MESH
%

disp([num2str(toc),'   GENERATING MESH'])

% MESH PROPERTIES

elemType  = 'Q4'; % the element type used in the FEM simulation;
elemTypeB = 'L2'; % the element type used in the FEM simulation;

%% domain 1 ---------------------------------------------------------------
numx1     = 10;
numy1     = 40;

L1 = 4.5;
% meshing for domain1
nnx=numx1+1;
nny=numy1+1;
%node1=square_node_array([0 0],[t 0],[t L1],[0 L1],nnx,nny);
node1=square_node_array([-c 0],[c 0],[c L1],[-c L1],nnx,nny);
%node1=square_node_array([-t 0],[0 0],[0 L1],[-t L1],nnx,nny);
inc_u=1;
inc_v=nnx;
node_pattern=[ 1 2 nnx+2 nnx+1 ];
element1 =make_elem(node_pattern,numx1,numy1,inc_u,inc_v);

%% domain 2 ---------------------------------------------------------------

numnode = 8;
node2 = zeros(numnode,2);
node2(:,1) = 0;
node2(:,2) = linspace(L1,L,numnode);

element2  = zeros(numnode-1,2);

for i=1:numnode-1
    element2(i,:) = i:i+1;
end

% PLOT MESH
if ( plotMesh )  
    clf
    plot_mesh(node1,element1,elemType,'g.-',2.4);
    plot_mesh(node2,element2,elemTypeB,'r.-',1.7);
    hold on
    axis off
    %axis([0 L -c c])
end

% boundary mesh for domain 1
bndMesh1 = numx1*numy1-numx1+1:numx1*numy1;

% boundary edges

bndNodes  = find(abs(node1(:,2)-L1)<1e-14);  
bndEdge1  = zeros(numx1,2);

for i=1:numx1
    bndEdge1(i,:) = bndNodes(i:i+1);
end


[W1,Q1]=quadrature( 2, 'GAUSS', 1 ); % two point quadrature

GP1 = [];
GP2 = [];
%% THIS NEEDS GENERALIZATION on y computation
for e=1:size(bndMesh1,2)
    sctrEdge = bndEdge1(e,:);
    sctr1    = element1(bndMesh1(e),:);    
    pts1     = node1(sctr1,:);    
    pts2     = node2(element2(1,:),:);   
    
    xx    = pts2(2,:) - pts2(1,:);
    cosphi = xx(1)/norm(xx);
    sinphi = xx(2)/norm(xx);
    Rv     = [cosphi sinphi;-sinphi cosphi];
    
    for q=1:size(W1)
        pt=Q1(q,:);
        wt=W1(q);
        [N,dNdxi]=lagrange_basis('L2',pt);  % element shape functions                
        J0    = dNdxi'*node1(sctrEdge,:);
        detJ0 = norm(J0);   
        J0    = J0/detJ0;        
        x     = N'*node1(sctrEdge,:); % global coord of GP
        xp  = pts2(1,:);
        yy  = Rv*(x-xp)';
        X1  = global2LocalMap(x,pts1,elemType);        
        GP1 = [GP1;X1 wt*detJ0 -J0(2) J0(1)];
        GP2 = [GP2;-1 yy(2)];
    end
end

% DEFINE BOUNDARIES

edgeElemType='L2';

% GET NODES ON DISPLACEMENT BOUNDARY
%      Here we get the nodes on the essential boundaries

fixedNode = find(node1(:,2)==0);  
midNode1  = find(node1(:,1)==0);  
rightNode = find(node2(:,2)==L);  

uFixed=zeros(1,length(fixedNode))';  % a vector of u_x for the nodes
vFixed=zeros(1,length(fixedNode))';     

numdofs = 2*size(node1,1) + 3*size(node2,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE SYSTEM DATA STRUCTURES

disp([num2str(toc),' INITIALIZING DATA STRUCTURES'])

f=zeros(numdofs,1);          % external load vector
K=zeros(numdofs,numdofs); % stiffness matrix

% ******************************************************************************
% ***                          P R O C E S S I N G                           ***
% ******************************************************************************
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%% COMPUTE STIFFNESS MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp([num2str(toc),'   COMPUTING STIFFNESS MATRIX'])

[W,Q]=quadrature( 2, 'GAUSS', 2 ); % 2x2 Gaussian quadrature

%% for domain 1------------------------------------------------------------

for e=1:size(element1,1)                          % start of element loop
    sctr = element1(e,:);           % element scatter vector
    nn   = length(sctr);
    sctrB(1:2:2*nn) = 2*sctr-1;
    sctrB(2:2:2*nn) = 2*sctr-0;    
    for q=1:size(W,1)
        pt=Q(q,:);
        wt=W(q);
        [N,dNdxi]=lagrange_basis(elemType,pt); % element shape functions
        J0=node1(sctr,:)'*dNdxi;                % element Jacobian matrix
        invJ0=inv(J0);
        dNdx=dNdxi*invJ0;
        % B matrix
        B(1,1:2:2*nn)  = dNdx(:,1)';
        B(2,2:2:2*nn)  = dNdx(:,2)';
        B(3,1:2:2*nn)  = dNdx(:,2)';
        B(3,2:2:2*nn)  = dNdx(:,1)';        
        % COMPUTE ELEMENT STIFFNESS AT QUADRATURE POINT
        K(sctrB,sctrB)=K(sctrB,sctrB)+B'*C*B*W(q)*det(J0);
    end  % of quadrature loop
end    

%% for domain 2 (beam)-----------------------------------------------------
clear sctrB;

[W,Q]=quadrature(  noGPs, 'GAUSS', 1 ); % 2 point quadrature

for e=1:size(element2,1)   
   conn  = element2(e,:);
   noFns = length(conn);   
   sctrB(1:3:3*noFns) = 3*conn-2;
   sctrB(2:3:3*noFns) = 3*conn-1;
   sctrB(3:3:3*noFns) = 3*conn-0;
   sctrB = sctrB + 2*size(node1,1);
   pts   = node2(conn,:);
   Ba    = zeros(1,noFns*3);
   
   xx    = pts(2,:) - pts(1,:);
   cosphi = xx(1)/norm(xx);
   sinphi = xx(2)/norm(xx);
   R(1:2,1:2) = [cosphi sinphi;-sinphi cosphi];
   R(4:5,4:5) = [cosphi sinphi;-sinphi cosphi];
   R(3,3) = 1;
   R(6,6) = 1;
   
   % loop over Gauss points 
    for gp=1:size(W,1)                        
      pt       = Q(gp,:);                          
      wt       = W(gp);                     
      [N,dNdxi]=lagrange_basis(elemTypeB,pt);  % element shape functions
      J0     = dNdxi'*pts;
      detJ0  = norm(J0);        
      dNdx   = (1/detJ0)*dNdxi;
                
      Ba(1:3:noFns*3)   = dNdx;
      Bb(3:3:noFns*3)   = dNdx;
      Bs(2:3:noFns*3)   = dNdx;
      Bs(3:3:noFns*3)   = -N;
      
      % compute elementary stiffness matrix       
      K(sctrB,sctrB) = K(sctrB,sctrB) + ...
          R'*(E0*A*Ba'*Ba + k * Bb' * Bb + shearFactor*G*A*Bs'*Bs)*R* detJ0 * wt;
    end
end

%% interface integrals-----------------------------------------------------

Cb     = [E0 0 0;
          0 0 0;
          0 0 shearFactor*G];
Cs = C;

Tinv  = zeros(3,3);
 
for i=1:numx1                     % start of element loop
    e1     = bndMesh1(i);    
    sctr1  = element1(e1,:);      
    sctr2  = element2(1,:); 
    %sctr2n = sctr2 + size(node1,1);
    nn1    = length(sctr1);
    nn2    = length(sctr2);
    
    sctrB1(1:2:2*nn1) = 2*sctr1-1; 
    sctrB1(2:2:2*nn1) = 2*sctr1-0; 
    
    sctrB2(1:3:3*nn2) = 3*sctr2-2; 
    sctrB2(2:3:3*nn2) = 3*sctr2-1; 
    sctrB2(3:3:3*nn2) = 3*sctr2-0; 
    
    sctrB2 = sctrB2 + 2*size(node1,1);
                
    pts1 = node1(sctr1,:);
    pts2 = node2(sctr2,:);
    
    Kp11 = zeros(nn1*2,nn1*2);
    Kp12 = zeros(nn1*2,nn2*3);
    Kp22 = zeros(nn2*3,nn2*3);
    
    Kd11 = zeros(nn1*2,nn1*2);
    Kd12 = zeros(nn1*2,nn2*3);
    Kd21 = zeros(nn2*3,nn1*2);
    Kd22 = zeros(nn2*3,nn2*3);
    
    xx    = pts2(2,:) - pts2(1,:);
    cosphi = xx(1)/norm(xx);
    sinphi = xx(2)/norm(xx);
    c2 = cosphi^2;
    s2 = sinphi^2;
    sc = sinphi*cosphi;
    
    R(1:2,1:2) = [cosphi sinphi;-sinphi cosphi];
    R(4:5,4:5) = [cosphi sinphi;-sinphi cosphi];
    R(3,3) = 1;
    R(6,6) = 1;
    Rv     = [cosphi sinphi;-sinphi cosphi];
    
    Tinv = [c2 s2 -2*sc;
            s2 c2 2*sc;
            sc -sc c2-s2];

    for q=2*i-1:2*i
        pt1 = GP1(q,1:2);
        wt1 = GP1(q,3);
        pt2 = GP2(q,1);
        y   = GP2(q,2);
        normal=GP1(q,4:5);
        n = [normal(1) 0 normal(2);
            0 normal(2) normal(1)];
        %n=-n;
        [N1,dN1dxi]=lagrange_basis(elemType,pt1); 
        [N2,dN2dxi]=lagrange_basis(elemTypeB,pt2); 
        
        J1 = pts1'*dN1dxi;                
        J2 = pts2'*dN2dxi;                
        detJ0=norm(J2);        
        
        dN1dx   = dN1dxi*inv(J1);        
        dN2dx   = (1/detJ0)*dN2dxi;
        
        Bs(1,1:2:2*nn1)  = dN1dx(:,1)';
        Bs(2,2:2:2*nn1)  = dN1dx(:,2)';
        Bs(3,1:2:2*nn1)  = dN1dx(:,2)';
        Bs(3,2:2:2*nn1)  = dN1dx(:,1)';
        
        Bb(1,1:3:3*nn2)  = dN2dx(:,1)';
        Bb(1,3:3:3*nn2)  = -y*dN2dx(:,1)';        
        Bb(3,2:3:3*nn2)  = dN2dx(:,1)';
        Bb(3,3:3:3*nn2)  = -N2;
        
        Ns(1,1:2:2*nn1)  = N1';
        Ns(2,2:2:2*nn1)  = N1';
        
        Nb(1,1:3:3*nn2)  = N2';
        Nb(1,3:3:3*nn2)  = -y*N2';
        Nb(2,2:3:3*nn2)  = N2';
        
        % COMPUTE ELEMENT STIFFNESS AT QUADRATURE POINT
        
        Kp11 = Kp11 + alpha*(Ns'*Ns)*wt1;
        Kp12 = Kp12 + alpha*(Ns'*Rv'*Nb*R)*wt1;
        Kp22 = Kp22 + alpha*(Nb'*Nb)*wt1;
        
        Kd11 = Kd11 + 0.5 * Ns'* n * Cs * Bs * wt1;
        Kd12 = Kd12 + 0.5 * Ns'* n * Tinv * Cb  * Bb * R *wt1;
        Kd21 = Kd21 + 0.5 * R'*Nb'* Rv * n * Cs * Bs *wt1;
        Kd22 = Kd22 + 0.5 * Nb'* Rv * n * Tinv*Cb  * Bb *wt1;
               
    end  % of quadrature loop
   
    Kd22 = R'*Kd22 * R;
    Kp22 = R'*Kp22 * R;
    
    K(sctrB1,sctrB1)  = K(sctrB1,sctrB1)  - Kd11 - Kd11' + Kp11;    
    K(sctrB1,sctrB2)  = K(sctrB1,sctrB2)  - Kd12 + Kd21' - Kp12;
    K(sctrB2,sctrB1)  = K(sctrB2,sctrB1)  + Kd21 - Kd12' - Kp12';
    K(sctrB2,sctrB2)  = K(sctrB2,sctrB2)  + Kd22 + Kd22' + Kp22;
end    


f(3*rightNode-2+ 2*size(node1,1)) = F;



%%%%%%%%%%%%%%%%%%% END OF STIFFNESS MATRIX COMPUTATION %%%%%%%%%%%%%%%%%%%


% APPLY ESSENTIAL BOUNDARY CONDITIONS
disp([num2str(toc),'   APPLYING BOUNDARY CONDITIONS'])
bcwt=mean(diag(K)); % a measure of the average size of an element in K
% used to keep the conditioning of the K matrix
udofs=2*fixedNode-1;           % global indecies of the fixed x displacements
vdofs=2*fixedNode;   % global indecies of the fixed y displacements
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

numnode1 = size(node1,1);
numnode2 = size(node2,1);

index1 = 1:numnode1;
index2 = 1:numnode2;

Ux1  = U(2*index1-1);
Uy1  = U(2*index1);

Ux2  = U(3*index2-2+2*numnode1);
Uy2  = U(3*index2-1+2*numnode1);
rot2 = U(3*index2-0+2*numnode1);

% Here we plot the stresses and displacements of the solution. As with the
% mesh generation section we don?t go into too much detail - use help
% ?function name? to get more details.
disp([num2str(toc),'   POST-PROCESSING'])

scaleFact=3.;
fn=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT DEFORMED DISPLACEMENT PLOT
%Ux = U(xs);
%Uy = U(ys);

colordef black
figure
clf
hold on
%plot_field(node1+scaleFact*[Ux1 Uy1],element1,elemType,Ux1);
%plot_field(node2+scaleFact*[Ux2 Uy2],element2,elemType,Ux2);
plot_mesh(node1+scaleFact*[Ux1 Uy1],element1,elemType,'w.-',1);
plot_mesh(node2+scaleFact*[Ux2 Uy2],element2,elemTypeB,'w.-',1);
%colorbar
axis off
%title('DEFORMED DISPLACEMENT IN Y-DIRECTION')

% Comapre numerical displacement to exact value
 
Uxm1  = U(2*midNode1-1);
Uym1  = U(2*midNode1);


xx = node1(midNode1,2);
u  = Uxm1;

xx=[xx;node2(:,2)];
u=[u;Ux2];

[W,Q]=quadrature(  10, 'GAUSS', 1 ); % 2 point quadrature

% for e=1:size(element2,1)   
%    conn  = element2(e,:);
%    noFns = length(conn);   
%    sctrx = 3*conn-2 + 2*size(node1,1);
%    
%    pts   = node2(conn,:);
%    
%    % loop over Gauss points 
%     for gp=1:size(W,1)                        
%       pt      = Q(gp,:);                          
%       wt      = W(gp);                     
%       [N,dNdxi]=lagrange_basis(elemTypeB,pt);  % element shape functions
%       J0=dNdxi'*pts;
%       xx = [xx; N'*pts(:,2)]; 
%       u  = [u;  N'*U(sctrx)]; 
%     end
% end

u = -u;

% exact solution u = 

x      = linspace(0,L,400);
uExact = -F/(6*k)*(3*L*x.^2-x.^3);

% solution from 2D continuum elements



colordef white
figure,set (gcf,'Color','w')
hold on
plot(x,uExact,'r-','LineWidth',2.);
plot(xx,u,'k-','MarkerSize',12,'LineWidth',1.4);
%plot(uCont(:,1),uCont(:,2),'b--','LineWidth',1.4);
%plot(uCoup(:,1),uCoup(:,2),'cy.-','LineWidth',1.4);
h=legend('exact','coupling','continuum','beam');
xlabel('x')
ylabel('w')
%axis([0 5 -0.55 0])

stress=zeros(size(element1,1),size(element1,2),3);

stressPoints=[-1 -1;1 -1;1 1;-1 1];


for e=1:size(element1,1)
    sctr=element1(e,:);
    sctrB(1:2:2*nn) = 2*sctr-1;
    sctrB(2:2:2*nn) = 2*sctr-0;
    nn=length(sctr);
    Ce=C;
  
    for q=1:nn
        pt=stressPoints(q,:);
        [N,dNdxi]=lagrange_basis(elemType,pt);
        J0=node1(sctr,:)'*dNdxi;
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


% sigmaXX = zeros(size(node1,1),2);
% sigmaYY = zeros(size(node1,1),2);
% sigmaXY = zeros(size(node1,1),2);
% 
% for e=1:size(element1,1)
%     connect = element1(e,:);
%     for in=1:4
%         nid = connect(in);
%         sigmaXX(nid,:) = sigmaXX(nid,:) + [stress(e,in,1) 1];
%         sigmaYY(nid,:) = sigmaYY(nid,:) + [stress(e,in,2) 1];
%         sigmaXY(nid,:) = sigmaXY(nid,:) + [stress(e,in,3) 1];               
%     end
% end
% 
% % Average nodal stress values (learned from Mathiew Pais XFEM code)
% sigmaXX(:,1) = sigmaXX(:,1)./sigmaXX(:,2); sigmaXX(:,2) = [];
% sigmaYY(:,1) = sigmaYY(:,1)./sigmaYY(:,2); sigmaYY(:,2) = [];
% sigmaXY(:,1) = sigmaXY(:,1)./sigmaXY(:,2); sigmaXY(:,2) = [];
% 
% plotStress(sigmaXX,sigmaXY,sigmaYY,element1,node1);

stressComp=2;
figure
clf
plot_field(node1+scaleFact*[Ux1 Uy1],element1,elemType,stress(:,:,stressComp));
plot_mesh(node2+scaleFact*[Ux2 Uy2],element2,elemTypeB,'b-',2);
%plot_mesh(node2+scaleFact*[Ux2 Uy2],element2,elemType,'r.-',1);
hold on
%plot_mesh(node+scaleFact*[U(xs) U(ys)],element,elemType,'g.-');
%plot_mesh(node,element,elemType,'w--');
colorbar

