% This file implements the Nitsche method to join two mechanical models:
% a 2D continuum model and a frame model.
% The beam elements can be placed arbitrarily in space.
% This is a prototype for solving frame problems with continuum elements
% in region of high gradients and beam elements in other places.
%
% Dofs: 
%     continuum nodes: ux and uy
%     beam nodes:      ux,uy and rotation
% Assembly: assume that continuum nodes are first followed by beam nodes.
%
% Discretisation: standard Lagrange finite elements.
%
% Problem: Timoshenko beam in bending.
% 
% Vinh Phu Nguyen
% Cardiff University, UK
% 7 June 2013

addpath ../fem_util/
addpath ../gmshFiles/
addpath ../post-processing/
addpath ../fem-functions/
addpath ../analytical-solutions/


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
L  = 5;        % length of the beam
c  = 0.5;      % the distance of the outer fiber of the beam from the mid-line

t  = 2*c;      % thickness, area = bxt
b  = 1; % 
I0=2*c^3/3;    % the second polar moment of inertia of the beam cross-section.
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

alpha = 600;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE FINITE ELEMENT MESH
%

disp([num2str(toc),'   GENERATING MESH'])

% MESH PROPERTIES

elemType  = 'Q4'; % the element type used in the FEM simulation;
elemTypeB = 'L2'; % the element type used in the FEM simulation;


numnode = 20;
node2 = zeros(numnode,2);
node2(:,1) = 0;
node2(:,2) = linspace(0,L,numnode);

element2  = zeros(numnode-1,2);

for i=1:numnode-1
    element2(i,:) = i:i+1;
end

% PLOT MESH
if ( plotMesh )  
    clf    
    plot_mesh(node2,element2,elemTypeB,'r.-',1.7);
    hold on
    axis off
    %axis([0 L -c c])
end

% DEFINE BOUNDARIES

edgeElemType='L2';

% GET NODES ON DISPLACEMENT BOUNDARY
%      Here we get the nodes on the essential boundaries

fixedNode=find(node2(:,2)==0);  
rightNode=find(node2(:,2)==L);  

uFixed=zeros(1,length(fixedNode))';  % a vector of u_x for the nodes
vFixed=zeros(1,length(fixedNode))';     
wFixed=zeros(1,length(fixedNode))';  

udofs = 3*fixedNode-2;
vdofs = 3*fixedNode-1;
wdofs = 3*fixedNode-0;

numdofs = 3*size(node2,1);

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



% for domain 2 (beam)
clear sctrB;

[W,Q]=quadrature(  noGPs, 'GAUSS', 1 ); % 2 point quadrature

for e=1:size(element2,1)   
   conn  = element2(e,:);
   noFns = length(conn);   
   sctrB(1:3:3*noFns) = 3*conn-2;
   sctrB(2:3:3*noFns) = 3*conn-1;
   sctrB(3:3:3*noFns) = 3*conn-0;   
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
          R'*(E0*A*Ba'*Ba + k * Bb' * Bb + shearFactor*G*b*t*Bs'*Bs)*R* detJ0 * wt;
    end
end

f(3*rightNode-2) = F;


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


numnode2 = size(node2,1);


index2 = 1:numnode2;



Ux2  = U(3*index2-2);
Uy2  = U(3*index2-1);
rot2 = U(3*index2-0);

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
plot_mesh(node2+scaleFact*[Ux2 Uy2],element2,elemTypeB,'w.-',1);
%colorbar
axis off
%title('DEFORMED DISPLACEMENT IN Y-DIRECTION')

% Comapre numerical displacement to exact value
 
xx = node2(:,2);
u  = -Ux2 ;


% exact solution u = 

x      = linspace(0,L,200);
uExact = -F/(6*k)*(3*L*x.^2-x.^3);

% solution from 2D continuum elements

uCont     = csvread('beam-2D.csv');
uCoup     = csvread('beam-1D.csv');

colordef white
figure,set (gcf,'Color','w')
hold on
plot(x,uExact,'r-','LineWidth',1.4);
plot(xx,u,'k-','MarkerSize',12,'LineWidth',1.4);
plot(uCont(:,1),uCont(:,2),'b--','LineWidth',1.4);
plot(uCoup(:,1),uCoup(:,2),'cy.-','LineWidth',1.4);
h=legend('exact','coupling','continuum','beam');
xlabel('x')
ylabel('w')
axis([0 5 -0.55 0])



