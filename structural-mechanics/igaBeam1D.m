%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Isogeometric analysis for one dimensional Euler-Bernoulli beam.
% (Rotation-free beam formulation due to high order continuity of 
% NURBS).
%
% Cantilever beam subjected to a point force.
%
% Vinh Phu Nguyen,
% Cardiff University
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

E = 100; % Young's modulus
t = 1; % thickness
b = 1;
l = 10; % length of the beam
I = b*t^3/12; % second moment of area
q = 0; % uniformly distributed loads
F = 1; % end force

k = E*I; 

% 1D beam [0,l]

refinement = 0;

knotVec    = [0 0 0 0 1 1 1 1];
controlPts = l*[0 0;0.3 0; 0.6 0; 1 0];
p          = 3;
noGPs      = 4;

knotVec = knotVec/max(knotVec);
weights = ones(1,size(controlPts,1));

% refinement and element connectivity 
generateIGA1DMesh


noCtrPts = size(controlPts,1);% no of control points
noElems  = size(elConn,1);    % no of elements

% initialization
K = zeros(noCtrPts,noCtrPts); % global stiffness matrix 
u = zeros(noCtrPts,1);        % displacement vector
f = zeros(noCtrPts,1);        % external force vector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gauss quadrature rule
[W,Q]=quadrature(  noGPs, 'GAUSS', 1 ); % 2 point quadrature

% Assembling system of equation
% Stiffness matrix and external force vector

disp([num2str(toc),'  ASSEMBLING THE SYSTEM'])

% Loop over elements (knot spans)
for e=1:noElems
   xiE   = elRange(e,:); % [xi_i,xi_i+1]
   conn  = elConn(e,:);
   noFns = length(conn);
   pts   = controlPts(conn,1:2);
   
   % loop over Gauss points 
    for gp=1:size(W,1)                        
      pt      = Q(gp,:);                          
      wt      = W(gp);                 
      % coord in parameter space  
      Xi      = 0.5 * (( xiE(2) - xiE(1) ) * pt + xiE(2) + xiE(1));
      J2      = 0.5 * ( xiE(2) - xiE(1) ); 
      
      % shape functions and first and second derivatives 
            
      [R dRdxi dR2dxi] = NURBS1DBasis2ndDers(Xi,p,knotVec,controlPts(:,3)');
      
      % compute the jacobian of physical and parameter domain mapping
      % then the derivative w.r.t spatial physical coordinates
      
      dxdxi   = dRdxi *pts;
      dx2dxi2 = dR2dxi*pts;
     
      J1     = norm (dxdxi);
      dRdx   = (1/J1)*dRdxi;
                
      dR2dx2 = 1/J1^2 * (dR2dxi - norm(dx2dxi2)*dRdx);
      
      % compute elementary stiffness matrix and
      % assemble it to the global matrix
      K(conn,conn) = K(conn,conn) + k * dR2dx2' * dR2dx2 * J1 * J2 * wt;
      
      % compute the external force, kind of body force
      f(conn) = f(conn) + q * R' * J1 * J2 * wt;
    end
end

%% Boundary conditions

% Dirichlet BCs
% fixed at the left end: no deflection, no rotation

udofs  = [1 2];  % global indecies  of the fixed x displacements
uFixed = [0 0]'; % BCs: u[0]=u[1]=0

% Neumman BCs
% concentrated force at the right end

f(noCtrPts) = F;

%% Solve the equation
bcwt=mean(diag(K)); % a measure of the average  size of an element in K
                    % used to keep the  conditioning of the K matrix
                                       
f=f-K(:,udofs)*uFixed;  % modify the  force vector

K(udofs,:)=0;  % zero out the rows and  columns of the K matrix
K(:,udofs)=0;
K(udofs,udofs)=bcwt*speye(length(udofs));  % put ones*bcwt on the diagonal
f(udofs)=bcwt*uFixed;

% SOLVE SYSTEM
disp([num2str(toc),'  SOLVING THE SYSTEM'])
U=K\f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% POST-PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
noXi = 20;
xi   = linspace(0,1,noXi);
pts  = zeros(noXi,1);
dis  = zeros(noXi,1);

for i=1:noXi
  dis(i,1) = NURBSinterpolation(xi(i), p, knotVec, ...
                     U', controlPts(:,3)');
  pts(i,1) = NURBSinterpolation(xi(i), p, knotVec, ...
                     controlPts(:,1)',controlPts(:,3)');
end

plot(pts,dis,'ks','MarkerSize',12);
hold on

% exact solution u = 

x      = linspace(0,l,120);
uExact = F/(6*k)*(3*l*x.^2-x.^3);

plot(x,uExact,'r-','LineWidth',1.4)
legend('IGA','exact')
set(gca,'Color',[0.8 0.8 0.8]);

% open the file with write permission










