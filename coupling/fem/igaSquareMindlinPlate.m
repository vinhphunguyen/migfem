%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FE analysis for Mindlin plate problems.
% Each node has three dofs(deflection,rota1,rota2).
%
% Fully clamped or simply supported rectangular plates.
%
% Nice example to illustrate shear locking for plate elements.
% Reduced integration as the easiest remedy: full integration for bending
% part and reduced integration for shear part.
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
addpath ../meshing/
addpath ../nurbs-util/
addpath ../nurbs-geopdes/inst/
addpath ../delamination/

clc
clear all

tic;

a = 1.0;
b = 1.0; 

% knots
uKnot = [0 0 1 1];
vKnot = [0 0 1 1];

% control points
controlPts          = zeros(4,2,2);

controlPts(1:2,1,1) = [0;0];
controlPts(1:2,2,1) = [a;0;];

controlPts(1:2,1,2) = [0;b];
controlPts(1:2,2,2) = [a;b];

% weights
controlPts(4,:,:)   = 1;

%% build NURBS object

solid = nrbmak(controlPts,{uKnot vKnot});

%% p-refinment

%solid = nrbdegelev(solid,[2 2]); % to cubic-cubic NURBS

%% h-refinement

refineLevel = 4;
for i=1:refineLevel
    uKnotVectorU = unique(uKnot);
    uKnotVectorV = unique(vKnot);
    
    % new knots along two directions
    
    newKnotsX = uKnotVectorU(1:end-1) + 0.5*diff(uKnotVectorU);
    newKnotsY = uKnotVectorV(1:end-1) + 0.5*diff(uKnotVectorV);
    
    newKnots  = {newKnotsX newKnotsY};
    solid     = nrbkntins(solid,newKnots);
    uKnot      = cell2mat(solid.knots(1));
    vKnot      = cell2mat(solid.knots(2));
end

mesh = buildIGA2DMesh (solid);

node     = mesh.controlPts;
element  = mesh.locElems;
element(:,[3 4]) = element(:,[4 3]);
elemType = 'Q4';

%% constitutive matrix

E  = 10920;
nu = 0.3;
t  = 0.1; % thickness

%% Boundary condition

q0  = -1.;  % distributed force

clamped = 0;

k   = 5/6;

D   = E*t^3/(12*(1-nu^2));
Cb  = D*[1 nu 0;nu 1 0; 0 0 0.5*(1-nu)];
Cs  = E*t*k/2/(1+nu)*[1 0;0 1];

% find boundary nodes for boundary conditions

EPS = 1e-8;
bottomNodes  =  find(abs(node(:,2))  <EPS);
topNodes     =  find(abs(node(:,2)-b)<EPS);
leftNodes    =  find(abs(node(:,1))  <EPS);
rightNodes   =  find(abs(node(:,1)-a)<EPS);

%fixedNodes   = leftNodes;
fixedNodes  =  unique([bottomNodes;rightNodes;topNodes;leftNodes]);

plot(node(fixedNodes,1),node(fixedNodes,2),...
    'bs','MarkerEdgeColor','r','MarkerSize',14);

% initialization

noElem = size(element,1);
noNode = size(node,1);
noDof  = 3*noNode;

K = sparse(noDof,noDof);  % global stiffness matrix
f = zeros(noDof,1);       % external force vector

% essential boundary conditions

uFixed     = zeros(size(fixedNodes))';
udofs      = 3*fixedNodes-2;

vFixed = [];
wFixed = [];
vdofs  = []; % rotation 1
wdofs  = []; % rotation 2 

if clamped
    vFixed  = zeros(size(fixedNodes))';
    wFixed  = zeros(size(fixedNodes))';
    vdofs   = 3*fixedNodes-1;
    wdofs   = 3*fixedNodes;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gauss quadrature rule
noGPs   = 2;
noGpEle = noGPs^2;
[W,Q]=quadrature(  noGPs, 'GAUSS', 2 ); % noGPs x noGPs point quadrature

% Assembling system of equation
% Stiffness matrix and external force vector

disp([num2str(toc),'  ASSEMBLING THE SYSTEM']);

% Loop over elements (knot spans)

for e=1:noElem
    sctr   = element(e,:);          %  element scatter vector
    nn     = length(sctr);
    pts    = node(sctr,:);
            
    nn3    = nn*3;
    sctrw  = 3*sctr-2;
    sctrB  = zeros(1,nn3);
    
    sctrB(1:3:nn3) = 3*sctr-2; % deflection
    sctrB(2:3:nn3) = 3*sctr-1; % rotation 1
    sctrB(3:3:nn3) = 3*sctr;   % rotation 2
    
    Bb     = zeros(3,3*nn);
    Bs     = zeros(2,3*nn);
    
    % loop over Gauss points
    
    for gp=1:size(W,1)
        pt      = Q(gp,:);
        wt      = W(gp);
        
        [R,dRdxi]= lagrange_basis(elemType,pt); % element shape functions
        J0       = pts'*dRdxi;                   % element Jacobian matrix
        invJ0    = inv(J0);
        dRdx     = dRdxi*invJ0;     
        
        % bending and shear B matrices
        
        Bb(1,2:3:3*nn) = dRdx(:,1);
        Bb(2,3:3:3*nn) = dRdx(:,2);
        Bb(3,2:3:3*nn) = dRdx(:,2);
        Bb(3,3:3:3*nn) = dRdx(:,1);
        
        Bs(1,2:3:3*nn) = -R;
        Bs(1,1:3:3*nn) = dRdx(:,1);
        Bs(2,1:3:3*nn) = dRdx(:,2);
        Bs(2,3:3:3*nn) = -R;
        
        % compute elementary stiffness matrix and
        % assemble it to the global matrix
        
        K(sctrB,sctrB) = K(sctrB,sctrB) + ...
                         (Bb' * Cb * Bb + Bs' * Cs * Bs) * det(J0)* wt;
        f(sctrw)       = f(sctrw) + q0 * R * det(J0) * wt;
    end
end

%f(3*rightNodes-2) = -0.5;

%%
disp([num2str(toc),'  APPLYING BOUNDARY  CONDITIONS']);

bcwt=mean(diag(K)); % a measure of the average  size of an element in K
% used to keep the  conditioning of the K matrix

% modify the force vector
f        = f-K(:,udofs)*uFixed';  
% f        = f-K(:,vdofs)*vFixed';  
% f        = f-K(:,wdofs)*wFixed';  
f(udofs) = bcwt*uFixed;
f(vdofs) = bcwt*vFixed;
f(wdofs) = bcwt*wFixed;

% modify the stiffness matrix
K(udofs,:)    = 0;  K(:,udofs)    = 0;
K(vdofs,:)    = 0;  K(:,vdofs)    = 0;
K(wdofs,:)    = 0;  K(:,wdofs)    = 0;

K(udofs,udofs)=bcwt*speye(length(udofs));  
K(vdofs,vdofs)=bcwt*speye(length(vdofs)); 
K(wdofs,wdofs)=bcwt*speye(length(wdofs)); 

%% SOLVE SYSTEM
disp([num2str(toc),'  SOLVING THE SYSTEM']);
U=K\f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% POST-PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ux      = U(2:3:noDof); % rotation1
Uy      = U(3:3:noDof); % rotation2
Uz      = U(1:3:noDof);

Ux(:) = 0;
Uy(:) = 0;

node3d = zeros(noNode,3);
node3d(:,[1 2]) = node;

scaleFact=1e0;
colordef black
figure
clf
hold on
plot_field(node3d+scaleFact*[Ux Uy Uz],element,elemType,Uz);
%plot_mesh(node3d+scaleFact*[Ux Uy Uz],element,elemType,'w.-',1);
colorbar
%axis off
%title('DEFORMED DISPLACEMENT IN Y-DIRECTION')
view(3)

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,fileName,opts)


% Exact solution check

wbar = min(Uz)*D/(a^4)
wext = 0.5  % CCCC
%wext = 4.06235; % SSSS







