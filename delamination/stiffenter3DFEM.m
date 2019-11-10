%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Isogeometric analysis for two dimensional elasticity problems.
% Using the so-called Bezier extraction operator.
%
% Cuvred composite panel with one stiffener. 3D version.
%
% Vinh Phu Nguyen,
% Cardiff University, UK
% April 2013.
% nvinhphu@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath ../fem_util/
addpath ../C_files/
addpath ../data/
addpath ../meshing/
addpath ../integration/
addpath ../post-processing/
addpath ../fem-functions/
addpath ../meshing/
addpath ../nurbs-util/
addpath ../nurbs-geopdes/inst/

clc

% input file
stiffener;


E0    = 3e7;  % Young modulus
nu0   = 0.25;  % Poisson ratio
g     = -3e4; % inner pressure

% COMPUTE ELASTICITY MATRIX
D=zeros(6,6);
D(1:3,1:3)=E0/(1+nu0)/(1-2*nu0)*[ 1-nu0 nu0 nu0;
    nu0 1-nu0 nu0;
    nu0 nu0 1-nu0];
D(4:6,4:6)=E0/2/(1+nu0)*eye(3);

tic;

% find boundary nodes for boundary conditions

fixedXNodes  =  data.xnodes'; % transpose to make it a row vector
fixedYNodes  =  data.ynodes';
fixedZNodes  =  data.znodes';

udofs = 3*fixedXNodes-2;    % global indecies  of the fixed x disps
vdofs = 3*fixedYNodes-1;    % global indecies  of the fixed y disps
wdofs = 3*fixedZNodes;      % global indecies  of the fixed z disps

uFixed = zeros(size(fixedXNodes));
vFixed = zeros(size(fixedYNodes));
wFixed = zeros(size(fixedZNodes));

% initialization
nDim       = 3;
patchCount = length(mesh);
noDofs     = data.pntCount * nDim;


K = sparse(noDofs,noDofs);  % global stiffness matrix
u = zeros(noDofs,1);        % displacement vector
f = zeros(noDofs,1);        % external force vector


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Assembling system of equation
% Stiffness matrix and external force vector

disp([num2str(toc),'  ASSEMBLING THE SYSTEM'])

% Loop over patches
for ip=1:patchCount
    mesh      = data.mesh{ip};
    globElems = mesh.globElems;
    locElems  = mesh.locElems;
    C         = mesh.C;
    weights   = mesh.weights;
    p         = mesh.p;
    q         = mesh.q;
    r         = mesh.r;
    controlPts= mesh.controlPts;
    
    % Gauss quadrature rule
    [W,Q] = gaussianQuadNURBS(p+1,q+1,r+1);
    
    %% Pre-compute Bernstein basis and derivatives for ONE Bezier element
    
    noBasis = (p+1)*(q+1)*(r+1);
    noGpEle = noBasis;
    
    shapes  = zeros(noGpEle,noBasis);
    derivs  = zeros(noGpEle,noBasis,3);
    
    for gp=1:size(W,1)
        [shapes(gp,:) derivs(gp,:,:)] = ...
            getShapeGradBernstein3D(p,q,r,Q(gp,1),Q(gp,2),Q(gp,3));
    end
    
    
    % Loop over elements
    for e=1:size(globElems,1)
        sctrg   = globElems(e,:);         %  global element scatter vector
        sctrl   = locElems (e,:);         %  local element scatter vector
        nn      = length(sctrg);
        
        sctrB(1:3:3*nn)    = 3*sctrg-2;
        sctrB(2:3:3*nn)    = 3*sctrg-1;
        sctrB(3:3:3*nn)    = 3*sctrg-0;
        B      = zeros(3,2*nn);
        
        Ce     = C(:,:,e);              % element Bezier extraction operator
        we     = diag(weights(sctrl));  % element weights
        pts    = controlPts(sctrl,:);   % element nodes
        Wb     = Ce'*weights(sctrl);    % element Bezier weights
        
        % loop over Gauss points
        for gp=1:size(W,1)
            pt      = Q(gp,:);
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
            
            if ip == 1
                sctry    = 3*sctrg-1;
                f(sctry) = f(sctry) + R * g * detJ * wt;
            end
        end
    end
end

w     = 1e10;
penaltyStiffness = w*[1 -1;-1 1];

for j=1:size(aa,2)
    for i=1:size(aa,1)
        sctr  = [aa(i,j) bb(i,j)];
        sctrx = 3*sctr-2;
        sctry = 3*sctr-1;
        sctrz = 3*sctr;
        K(sctrx,sctrx) = K(sctrx,sctrx) + penaltyStiffness;
        K(sctry,sctry) = K(sctry,sctry) + penaltyStiffness;
        K(sctrz,sctrz) = K(sctrz,sctrz) + penaltyStiffness;
    end
end

%% Computing external force

disp([num2str(toc),'  APPLYING BOUNDARY  CONDITIONS'])

[K,f]  = applyDirichletBCs3D(K,f,udofs,vdofs,wdofs,uFixed,vFixed,wFixed);

% SOLVE SYSTEM
disp([num2str(toc),'  SOLVING THE SYSTEM'])
U=K\f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% POST-PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


vtuFile0 = '../results/stiffener';


% Loop over patches
for ip=1:patchCount    
    vtuFile = strcat(vtuFile0,num2str(ip));
    figure; hold on;
    ok      = plotStress3DForPatch(data,ip,vtuFile,U,D);
end










