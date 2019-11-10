%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Isogeometric analysis for one dimensional Timoshenko beam.
% Each control point has two dofs (deflection and rotation).
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

E  = 1000; % Young's modulus
nu = 0.3; % Poisson ratio
t  = 1; % thickness, area = bxt
b  = 1; % 
l  = 5; % length of the beam
I  = b*t^3/12; % second moment of area
q  = 0; % uniformly distributed loads
F  = -1; % end force
G  = E/2/(1+nu);
k  = E*I; 
shearFactor = 5/6;

% 1D beam [0,l]

refinement = 2;

knotVec    = [0 0 0 0 1 1 1 1];
controlPts = l*[0 0;0.3 0; 0.6 0; 1 0];
p          = 3;
noGPs      = 4;

knotVec = knotVec/max(knotVec);
weights = ones(1,size(controlPts,1));

% refinement and element connectivity 
generateIGA1DMesh


noCtrPts = size(controlPts,1);% no of control points
noDofs   = noCtrPts * 2;
noElems  = size(elConn,1);    % no of elements

% initialization
K = zeros(noDofs,noDofs); % global stiffness matrix 
u = zeros(noDofs,1);        % displacement vector
f = zeros(noDofs,1);        % external force vector

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
   sctrB(1:2:2*noFns) = 2*conn-1;
   sctrB(2:2:2*noFns) = 2*conn-0;
   
   pts   = controlPts(conn,1:2);
   
   % loop over Gauss points 
    for gp=1:size(W,1)                        
      pt      = Q(gp,:);                          
      wt      = W(gp);                 
      % coord in parameter space  
      Xi      = 0.5 * (( xiE(2) - xiE(1) ) * pt + xiE(2) + xiE(1));
      J2      = 0.5 * ( xiE(2) - xiE(1) ); 
      
      % shape functions and first derivatives 
            
      [R dRdxi] = NURBS1DBasisDers(Xi,p,knotVec,controlPts(:,3)');
      
      % compute the jacobian of physical and parameter domain mapping
      % then the derivative w.r.t spatial physical coordinates
      
      dxdxi   = dRdxi *pts;      
     
      J1     = norm (dxdxi);
      dRdx   = (1/J1)*dRdxi;
                
      Bb(2:2:noFns*2)   = dRdx;
      Bs(1:2:noFns*2)   = dRdx;
      Bs(2:2:noFns*2)   = -R;
      
      % compute elementary stiffness matrix and
      % assemble it to the global matrix
      K(sctrB,sctrB) = K(sctrB,sctrB) + ...
          (k * Bb' * Bb + shearFactor*G*b*t*Bs'*Bs)* J1 * J2 * wt;
      
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

f(2*noCtrPts-1) = F;

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
                     U(1:2:noDofs)', controlPts(:,3)');
  pts(i,1) = NURBSinterpolation(xi(i), p, knotVec, ...
                     controlPts(:,1)',controlPts(:,3)');
end

% exact solution u = 

x      = linspace(0,l,120);
uExact = F/(6*k)*(3*l*x.^2-x.^3);

% solution from 2D continuum elements

data = [pts dis];
csvwrite('beam-1D.csv',data);










