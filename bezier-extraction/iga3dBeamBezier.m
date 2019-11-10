%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Isogeometric analysis for three dimensional elasticity problems.
%
% 3D beam in bending.
% Illustration of 3D Bezier extraction.
%
% Vinh Phu Nguyen,
% Cardiff University, Wales, UK
% nvinhphu@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath ../fem_util/
addpath ../C_files/
addpath ../data/
addpath ../meshing/
addpath ../post-processing/
addpath ../fem-functions/
addpath ../nurbs-util/
addpath ../integration/
addpath ../analytical-solutions/

clc
clear all

global p q r

E0           = 1e5;  % Young modulus
nu0          = 0.3;  % Poisson’s ratio

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE COMPLIANCE MATRIX
D=zeros(6,6);
D(1:3,1:3)=E0/(1+nu0)/(1-2*nu0)*[ 1-nu0 nu0 nu0;
    nu0 1-nu0 nu0;
    nu0 nu0 1-nu0];
D(4:6,4:6)=E0/2/(1+nu0)*eye(3);

tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%  CAD input: control points and knot vectors
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


beam3dC1Data
noGPs        = 3; % # of Gauss points along one direction

[C,Cxi,Cet,Cze] = bezierExtraction3D(uKnot,vKnot,wKnot,p,q,r);

noCtrPts   = noPtsX * noPtsY * noPtsZ;
noDofs     = noCtrPts * 3;

% Boundary nodes
% 1. Displacement conditions:
%    - On bottom edge (y=0) with u_y = 0
%    - On right  edge (x=0) with u_x = 0

% 2. Traction conditions
%    - On left edge (x=-L) with t = (sigma_x,sigma_xy)
%    - On top edge  (y=L) with t = (sigma_xy,sigma_y)

% find boundary nodes for bounjdary conditions

leftNodes   = find(controlPts(:,1)==0)';
rightNodes  = find(controlPts(:,1)==a)';

% essential boundary conditions

uFixed     = zeros(size(leftNodes));
vFixed     = zeros(size(leftNodes));
wFixed     = zeros(size(leftNodes));

% build connectivity ...

generateIGA3DMesh

% initialization

K = sparse(noDofs,noDofs);  % global stiffness matrix
u = zeros(noDofs,1);        % displacement vector
f = zeros(noDofs,1);        % external force vector

%% Pre-compute Bernstein basis and derivatives for ONE Bezier element

[W,GP] = gaussianQuadNURBS(p+1,q+1,r+1);

noBasis = (p+1)*(q+1)*(r+1);
noGpEle = noGPs*noGPs*noGPs;

shapes  = zeros(noGpEle,noBasis);
derivs  = zeros(noGpEle,noBasis,3);

for gp=1:size(W,1)
    [shapes(gp,:) derivs(gp,:,:)] = ...
        getShapeGradBernstein3D(p,q,r,GP(gp,1),GP(gp,2),GP(gp,3));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Assembling system of equation
% Stiffness matrix and external force vector

disp([num2str(toc),'  ASSEMBLING THE SYSTEM'])

% Loop over elements (knot spans)

for e=1:noElems
    idu    = index(e,1);
    idv    = index(e,2);
    idw    = index(e,3);
    
    xiE    = elRangeU(idu,:); % [xi_i,xi_i+1]
    etaE   = elRangeV(idv,:); % [eta_j,eta_j+1]
    zetaE  = elRangeW(idw,:); % [zeta_k,zeta_k+1]
    
    sctr   = element(e,:);          %  element scatter vector    
    nn     = length(sctr);
    
    sctrB(1:3:3*nn)    = 3*sctr-2;
    sctrB(2:3:3*nn)    = 3*sctr-1;
    sctrB(3:3:3*nn)    = 3*sctr-0;
    
    Ce     = C(:,:,e);             % element Bezier extraction operator
    we     = diag(weights(sctr));  % element weights
    pts    = controlPts(sctr,:);   % element nodes
    Wb     = Ce'*weights(sctr);    % element Bezier weights
        
    % loop over Gauss points
    
    for gp=1:size(W,1)        
        wt      = W(gp);
        
        %% Bernstein basis and derivatives at GP gp
        Be      = shapes(gp,:)';
        dBedxi  = reshape(derivs(gp,:,:),noBasis,3);
        
        %% Bezier weight functions (denomenator of NURBS)
        wb        = dot(Be,Wb);            % Be(I)*Wb(I)
        dwbdxi(1) = dot(dBedxi(:,1),Wb);   % Be(I)_{,xi} * Wb(I)
        dwbdxi(2) = dot(dBedxi(:,2),Wb);   % Be(I)_{,et} * Wb(I)
        dwbdxi(3) = dot(dBedxi(:,3),Wb);   % Be(I)_{,et} * Wb(I)
        %% Shape function and derivatives
        R          = we*Ce*Be/wb;
        dRdxi(:,1) = we*Ce*(dBedxi(:,1)/wb-dwbdxi(1)*Be/(wb*wb));
        dRdxi(:,2) = we*Ce*(dBedxi(:,2)/wb-dwbdxi(2)*Be/(wb*wb));
        dRdxi(:,3) = we*Ce*(dBedxi(:,3)/wb-dwbdxi(3)*Be/(wb*wb));
           
        %% Jacobian matrix
        dxdxi = pts'*dRdxi;
        
        dxidx = inv(dxdxi);
        dRdx  = dRdxi*dxidx;
        detJ  = det(dxdxi);
        
        % B matrix
        
        B = strainDispMatrix3d(nn,dRdx);
        
        % compute elementary stiffness matrix and
        % assemble it to the global matrix
        
        K(sctrB,sctrB) = K(sctrB,sctrB) + B' * D * B * detJ * wt;
    end
end

% Computing external force

f(3*rightNodes)=-100;

disp([num2str(toc),'  APPLYING BOUNDARY  CONDITIONS'])

udofs = 3*leftNodes-2;    % global indecies  of the fixed x disps
vdofs = 3*leftNodes-1;    % global indecies  of the fixed y disps
wdofs = 3*leftNodes;      % global indecies  of the fixed z disps

[K,f]  = applyDirichletBCs3D(K,f,udofs,vdofs,wdofs,uFixed,vFixed,wFixed);

% SOLVE SYSTEM

disp([num2str(toc),'  SOLVING THE SYSTEM'])
U=K\f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% POST-PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp([num2str(toc),'  POST-PROCESSING'])

Ux    = U(1:3:noDofs);
Uy    = U(2:3:noDofs);
Uz    = U(3:3:noDofs);

vtsFile = '../results/beam3D';
plotStress3d

figure
fac=2;
view(3)
hold on
plot_mesh(node+fac*[dispX dispY dispZ],elementV,'B8','r.-',1.2);


% figure
% hold on
% plot_mesh(node,elementV,'B8','g.-');
% %plot_field(node,elementV,'B8',Uz); do not support B8 elements!!!
% view(3)
% plot3(controlPts(:,1),controlPts(:,2),controlPts(:,3),'r*');
%
% buildVisualization3dMesh
% sigmaXX = zeros(size(node,1),1);
% sigmaYY = zeros(size(node,1),1);
% sigmaXY = zeros(size(node,1),1);
%
% VTKPostProcess3d(node,elementV,'B8',...
%              [sigmaXX sigmaYY sigmaXY sigmaXX sigmaXX sigmaXX],[Ux Uy Uz]);



