%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Finite elements for two dimensional frame analysis.
% Frame element = beam + bar, each node has three unknowns:
% two displacements and one rotation. 
% Element stiffness in local coordinate system and transformed to global
% one using the so-called rotation matrix.
%
% Vinh Phu Nguyen,
% Cardiff University, 7 June 2013
% nvinhphu@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath ../fem_util/
addpath ../C_files/
addpath ../data/
addpath ../meshing/
addpath ../post-processing/
addpath ../fem-functions/
addpath ../nurbs-util/

clc
clear all

global p

tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% PRE-PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Material, geometry and force data

E  = 1000; % Young's modulus
nu = 0.3; % Poisson ratio

b = 1;
t = 1;

q  = -1; % uniformly distributed loads
G  = E/2/(1+nu);
A  = b*t;
I  = b*t^3/12;
k  = E*I; 
shearFactor = 5/6;

%% Mesh
l = 10;

node    = [0 0;0 l;l/2 l; l l;l 0];
element = [1 2;2 3;3 4;4 5];
elemType = 'L2';

% force
F = 1;
force = [0 0 0 0];

noCtrPts = size(node,1);% no of control points
noDofs   = noCtrPts * 3;
noElems  = size(element,1);    % no of elements

%% initialization
K = zeros(noDofs,noDofs); % global stiffness matrix 
u = zeros(noDofs,1);        % displacement vector
f = zeros(noDofs,1);        % external force vector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gauss quadrature rule
[W,Q]=quadrature(  2, 'GAUSS', 1 ); % 2 point quadrature

% Assembling system of equation
% Stiffness matrix and external force vector

disp([num2str(toc),'  ASSEMBLING THE SYSTEM'])

% Loop over elements (knot spans)
R = zeros(6,6);
for e=1:noElems   
   conn  = element(e,:);
   noFns = length(conn);
   sctrB(1:3:3*noFns) = 3*conn-2;
   sctrB(2:3:3*noFns) = 3*conn-1;
   sctrB(3:3:3*noFns) = 3*conn-0;
   sctrf = 3*conn-1; % assume force in the y direction
   pts   = node(conn,:);
   q     = force(e);
   xx    = pts(2,:) - pts(1,:);
   cosphi = xx(1)/norm(xx);
   sinphi = xx(2)/norm(xx);
   R(1:2,1:2) = [cosphi sinphi;-sinphi cosphi];
   R(4:5,4:5) = [cosphi sinphi;-sinphi cosphi];
   R(3,3) = 1;
   R(6,6) = 1;
   % rotation matrix
   Ba = zeros(1,noFns*3);
   % loop over Gauss points 
    for gp=1:size(W,1)                        
      pt      = Q(gp,:);                          
      wt      = W(gp);                         
      % shape functions and first derivatives             
      [N,dNdxi]=lagrange_basis(elemType,pt);  % element shape functions
      J0=dNdxi'*pts;
      detJ0=norm(J0);
      dNdx   = (1/detJ0)*dNdxi;
      
      Ba(1:3:noFns*3)   = dNdx;
      Bb(3:3:noFns*3)   = dNdx;
      Bs(2:3:noFns*3)   = dNdx;
      Bs(3:3:noFns*3)   = -N;
      
      % compute elementary stiffness matrix       
      K(sctrB,sctrB) = K(sctrB,sctrB) + ...
          R'*(E*A*Ba'*Ba + k * Bb' * Bb + shearFactor*G*b*t*Bs'*Bs)* R* detJ0 * wt;
      
      % compute the external force, kind of body force
      f(sctrf) = f(sctrf) + q * N * detJ0* wt;
    end
end

f(3*3-1) = -F;

%% Boundary conditions

fixedNode = [1 5];

udofs     = 3*fixedNode-2;
vdofs     = 3*fixedNode-1;
wdofs     = 3*fixedNode-0;

uFixed    = zeros(1,length(fixedNode));
vFixed    = zeros(1,length(fixedNode));
wFixed    = zeros(1,length(fixedNode));

% Neumman BCs
% concentrated force at the right end


%% Solve the equation

[K,f]=applyDirichletBCs3D(K,f,udofs,vdofs,wdofs,uFixed,vFixed,wFixed);

% SOLVE SYSTEM
disp([num2str(toc),'  SOLVING THE SYSTEM'])
U=K\f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% POST-PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

Ux      = U(1:3:noDofs);
Uy      = U(2:3:noDofs);


% Here we plot the stresses and displacements of the solution.
disp([num2str(toc),'   POST-PROCESSING'])

dispNorm=l/max(sqrt(Ux.^2+Uy.^2));
scaleFact=0.1*dispNorm;

colordef white
figure
clf
hold on
%plot_field(node1+scaleFact*[Ux1 Uy1 Uz1],element1,elemType,Uy1);
%plot_field(node2+scaleFact*[Ux2 Uy2 Uz2],element2,elemType,Uy2);
plot_mesh(node+scaleFact*[Ux Uy],element,elemType,'r-',1);
%colorbar
axis off






