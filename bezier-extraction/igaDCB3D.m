%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Isogeometric analysis for three dimensional delamination problems.
%
% 3D double cantilever beam.
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
f0           = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE COMPLIANCE MATRIX
D=zeros(6,6);
D(1:3,1:3)=E0/(1+nu0)/(1-2*nu0)*[ 1-nu0 nu0 nu0;
    nu0 1-nu0 nu0;
    nu0 nu0 1-nu0];
D(4:6,4:6)=E0/(1+nu0)*eye(3);

tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%  CAD input: control points and knot vectors
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


dcbCk3DData
noGPs        = 3; % # of Gauss points along one direction

noCtrPts   = noPtsX * noPtsY * noPtsZ;
noDofs     = noCtrPts * 3;


% find boundary nodes for bounjdary conditions

eps = 1e-10;

rightNodes  = find(abs(controlPts(:,1)-L) < eps)';
leftNodes   = find(abs(controlPts(:,1)-0) < eps)';
botNodes    = find(abs(controlPts(:,2)+W) < eps)';
topNodes    = find(abs(controlPts(:,2)-W) < eps)';

forcedNodes1 = intersect(leftNodes,botNodes);
forcedNodes2 = intersect(leftNodes,topNodes);

% essential boundary conditions

uFixed     = zeros(size(rightNodes));
vFixed     = zeros(size(rightNodes));
wFixed     = zeros(size(rightNodes));

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
    sctrB  = [sctr sctr+noCtrPts sctr+2*noCtrPts]; % scatters a B matrix
    nn     = length(sctr);
    
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

f(forcedNodes1+noCtrPts) = -f0;
f(forcedNodes2+noCtrPts) =  f0;

disp([num2str(toc),'  APPLYING BOUNDARY  CONDITIONS'])

udofs = rightNodes;             % global indecies  of the fixed x disps
vdofs = rightNodes+noCtrPts;    % global indecies  of the fixed y disps
wdofs = rightNodes+2*noCtrPts;  % global indecies  of the fixed z disps

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

Ux    = U(1:noCtrPts);
Uy    = U(noCtrPts+1:2*noCtrPts);
Uz    = U(2*noCtrPts+1:noDofs);

%% convert to NURBS object for visualization

fac = 1e3;
defControlPts = controlPts+fac*[Ux Uy Uz];

controlPts1          = zeros(4,noPtsX,noPtsY,noPtsZ);

noPtsXY = noPtsX * noPtsY;

for i=1:noPtsZ
    for j=1:noPtsY
        en  = noPtsXY*(i-1) + noPtsX*(j-1) + 1;
        st  = en + noPtsX - 1;
        controlPts1(1:3,1:noPtsX,j,i) = defControlPts(en:st,:)';
    end
end

controlPts1(4,:,:)   = 1;

solid = nrbmak(controlPts1,{uKnot vKnot wKnot});

nrbplot(solid,[20 20 4]);




