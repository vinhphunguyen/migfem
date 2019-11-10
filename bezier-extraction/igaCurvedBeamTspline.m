%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Isogeometric analysis for two dimensional elasticity problems.
% Tspline geometry.
% Using the so-called Bezier extraction operator.
% Using the triple sparse matrix format to speed up the assembly process.
%
% Cylinder subject to inner pressure. Only a quarter is modeled.
%
% Vinh Phu Nguyen,
% Cardiff University, UK
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

clc
clear all

global p q

E0          = 3e7;           % Young modulus
nu0         = 0.25;          % Poisson ratio
stressState ='PLANE_STRESS'; % either 'PLANE_STRAIN' or "PLANE_STRESS
pressure    = 3e4;           % inner pressure

% COMPUTE ELASTICITY MATRIX

De = elasticityMatrix(E0,nu0,stressState);

tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%  CAD input: *.iga file from Rhino3d plugin
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filename = 'arch_tsplines.iga';
tspline  = read_bezier_extraction (filename);
convert2DTsplines;

% find boundary nodes for Dirichlet boundary conditions
EPS = 1e-8;
fixedXNodes  =  find(abs(controlPts(:,1))<EPS);
fixedYNodes  =  find(abs(controlPts(:,2))<EPS);

%% initialization

u = zeros(noDofs,1);        % displacement vector
f = zeros(noDofs,1);        % external force vector

%% one the assembly
nElNod = nsh_max;
nElDof = nElNod*2;
nElmLK = nElDof^2;
nSprGK = nElmLK*noElems;

jSprLK = 1:nElmLK;
vIdGrd = ones(1,nElDof);

% values, row indices, columns indices of the global K matrix
vSprGK = zeros(nSprGK,1);
jSprRw = zeros(nSprGK,1);
jSprCl = zeros(nSprGK,1);

Ke0    = zeros(nElDof,nElDof); % element Ke

%% essential boundary conditions

uFixed     = zeros(size(fixedXNodes))';
vFixed     = zeros(size(fixedYNodes))';

udofs      = 2*fixedXNodes-1;
vdofs      = 2*fixedYNodes;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gauss quadrature rule
noGPs       = max(p,q)+1; % # of Gauss points along one direction
[W,Q]=quadrature(  noGPs, 'GAUSS', 2 ); % noGPs x noGPs point quadrature

%% Pre-compute Bernstein basis and derivatives for ONE Bezier element

noBasis = (p+1)*(q+1);
noGpEle = noGPs*noGPs;

shapes  = zeros(noGpEle,noBasis);
derivs  = zeros(noGpEle,noBasis,2);

for gp=1:size(W,1)
    [shapes(gp,:) derivs(gp,:,:)] = getShapeGradBernstein2D(p,q,Q(gp,1),Q(gp,2));
end

%% Assembling system of equation
% Stiffness matrix and external force vector

disp([num2str(toc),'  ASSEMBLING THE SYSTEM'])

% Loop over elements

for e=1:noElems
    sctr   = element{e};           %  element scatter vector    
    nn     = length(sctr);
    
    sctrB(1,1:2:2*nn) = 2*sctr-1;
    sctrB(1,2:2:2*nn) = 2*sctr;
    
    B      = zeros(3,2*nn);
    Ke     = Ke0;
    
    Ce     = C{e};                 % element Bezier extraction operator
    we     = diag(weights(sctr));  % element weights
    pts    = controlPts(sctr,:);   % element nodes
    Wb     = Ce'*weights(sctr);    % element Bezier weights
    
    % loop over Gauss points
    
    for gp=1:size(W,1)
        %pt      = Q(gp,:);
        wt      = W(gp);
        
        %% Bernstein basis and derivatives at GP gp
        Be      = shapes(gp,:)';
        dBedxi  = reshape(derivs(gp,:,:),noBasis,2);
        
        %% Bezier weight functions (denomenator of NURBS)
        wb        = dot(Be,Wb);
        dwbdxi(1) = dot(dBedxi(:,1),Wb);
        dwbdxi(2) = dot(dBedxi(:,2),Wb);
        %% Shape function and derivatives
        R          = we*Ce*Be/wb;
        dRdxi(:,1) = we*Ce*(dBedxi(:,1)/wb-dwbdxi(1)*Be/(wb*wb));
        dRdxi(:,2) = we*Ce*(dBedxi(:,2)/wb-dwbdxi(2)*Be/(wb*wb));
        
        %% Jacobian matrix
        dxdxi = pts'*dRdxi;
        
        dxidx = inv(dxdxi);
        dRdx  = dRdxi*dxidx;
        detJ  = det(dxdxi);
        
        % B matrix
        B(1,1:2:2*nn)  = dRdx(:,1)';
        B(2,2:2:2*nn)  = dRdx(:,2)';
        B(3,1:2:2*nn)  = dRdx(:,2)';
        B(3,2:2:2*nn)  = dRdx(:,1)';
        
        % compute elementary stiffness matrix and
        % assemble it to the global matrix
        
        Ke = Ke + B' * De * B * detJ * wt;
    end
    
    %----------------------------------------------------------------------
    % Global Assembly (Sparse)
    %----------------------------------------------------------------------
    
    vElDof = sctrB';
    mRwGrd = vElDof(:,vIdGrd);
    mClGrd = mRwGrd';
    
    
    jSprRw(jSprLK) = mRwGrd(:);
    jSprCl(jSprLK) = mClGrd(:);
    vSprGK(jSprLK) = Ke(:);
    
    jSprLK         = jSprLK + nElmLK; % move to the next position
end


% Here comes the total stiffness matrix in one shot!!!

K = sparse(jSprRw,jSprCl,vSprGK,noDofs,noDofs);

%% Computing external force

[W1,Q1] = quadrature(noGPs, 'GAUSS', 1 );


%% find edges for Neumann BCs

ibnd     = 1; % inner surface
bndElems = tspline.elements(tspline.boundary{ibnd}.elem);
nsides   = tspline.boundary{ibnd}.nsides;

for iside = 1:nsides
    bndElem = bndElems(iside);
    degree  = bndElem.degree;
    % Identify the univariate Bernstein polynomial from the position,
    %   L(eft), R(ight), B(ottom) and T(op).
    switch tspline.boundary{ibnd}.position{iside}
        case {'L'}
            bernstein_indices = sub2ind (degree+1, ones (1, degree(2)+1), 1:degree(2)+1);
            ind = 2;
        case {'R'}
            bernstein_indices = sub2ind (degree+1, (degree(1)+1)*ones (1, degree(2)+1), 1:degree(2)+1);
            ind = 2;
        case {'B'}
            bernstein_indices = sub2ind (degree+1, 1:degree(1)+1, ones (1, degree(1)+1));
            ind = 1;
        case {'T'}
            bernstein_indices = sub2ind (degree+1, 1:degree(1)+1, (degree(2)+1)*ones (1, degree(1)+1));
            ind = 1;
    end
    
    % Compute the reduced extraction operator
    bndC        = bndElem.extraction;
    bsp_indices = any (bndC(:,bernstein_indices), 2);
    C_bnd       = bndC(bsp_indices, bernstein_indices);%1D Bezier extraction
    sctr        = bndElem.connectivity(bsp_indices);
    sctrx       = 2*sctr-1;
    sctry       = 2*sctr;
    nsh_iside   = numel (sctrx);
    
    pts         = controlPts(sctr,:);
    wes         = weights(sctr);        
    we          = diag(wes);               % element weights
    Wb          = C_bnd'*wes;               % element Bezier weights
        
    % loop over Gauss points
    for gp=1:size(W1,1)
        xi      = Q1(gp,:);
        wt      = W1(gp);
        
        [Be,dBdxi]= bernsteinBasis1D(q,xi);
        wb        = dot(Be,Wb);
        R         = we*C_bnd*Be/wb;
        dwbdxi    = dot(dBdxi,Wb);
        dRdxi     = we*C_bnd*(dBdxi/wb-dwbdxi*Be/(wb*wb));
        
        x        = R' * pts; % global coord of GP
        
        jacob    = dRdxi' * pts;
        J        = norm (jacob);
        
        r        = norm(x);
        Fx       = pressure * x(1,1)/r;
        Fy       = pressure * x(1,2)/r;
        
        f(sctrx) = f(sctrx) + R * Fx * J * wt;
        f(sctry) = f(sctry) + R * Fy * J * wt;
    end
end
dRdxi=[];
%%
disp([num2str(toc),'  APPLYING BOUNDARY  CONDITIONS'])

[K,f]=applyDirichletBCs(K,f,udofs,vdofs,uFixed,vFixed);


% SOLVE SYSTEM
disp([num2str(toc),'  SOLVING THE SYSTEM'])
U=K\f;


%%
%%%
%%% POST-PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ux    = U(1:2:noDofs);
Uy    = U(2:2:noDofs);

dRdxi=[];
x     = zeros(3,noGpEle,noElems);  % global coords of Gauss points
u     = zeros(3,noGpEle,noElems);  % displacements of Gauss points
sigma = zeros(3,noGpEle,noElems);  % stresses      of Gauss points

id = 1;
xx=zeros(noGpEle*noElems,2);
for e=1:noElems
    sctr   = element{e};           %  element scatter vector    
    nn     = length(sctr);
    
    sctrB(1,1:2:2*nn) = 2*sctr-1;
    sctrB(1,2:2:2*nn) = 2*sctr;
    
    B      = zeros(3,2*nn);
    Ke     = Ke0;
    
    Ce     = C{e};                 % element Bezier extraction operator
    we     = diag(weights(sctr));  % element weights
    pts    = controlPts(sctr,:);   % element nodes
    Wb     = Ce'*weights(sctr);    % element Bezier weights
    
    elemDisp = U(sctrB);
    
    % loop over Gauss points
    
    for gp=1:size(W,1)
        pt      = Q(gp,:);
        wt      = W(gp);
        
        %% Bernstein basis and derivatives at GP gp
        Be      = shapes(gp,:)';
        dBedxi  = reshape(derivs(gp,:,:),noBasis,2);
        
        %% Bezier weight functions (denomenator of NURBS)
        wb        = dot(Be,Wb);
        dwbdxi(1) = dot(dBedxi(:,1),Wb);
        dwbdxi(2) = dot(dBedxi(:,2),Wb);
        %% Shape function and derivatives
        R          = we*Ce*Be/wb;
        dRdxi(:,1) = we*Ce*(dBedxi(:,1)/wb-dwbdxi(1)*Be/(wb*wb));
        dRdxi(:,2) = we*Ce*(dBedxi(:,2)/wb-dwbdxi(2)*Be/(wb*wb));
        
        %% Jacobian matrix
        dxdxi = pts'*dRdxi;
        
        dxidx = inv(dxdxi);
        dRdx  = dRdxi*dxidx;
        detJ  = det(dxdxi);
        
        % B matrix
        B(1,1:2:2*nn)  = dRdx(:,1)';
        B(2,2:2:2*nn)  = dRdx(:,2)';
        B(3,1:2:2*nn)  = dRdx(:,2)';
        B(3,2:2:2*nn)  = dRdx(:,1)';
        
        x(1:2,gp,e)     = R' * pts;
        xx(id,:)        = R' * pts;
        u(1:2,gp,e)     = R' * [Ux(sctr) Uy(sctr)];
        sigma(:,gp,e)   = De * B * elemDisp;
        id = id + 1;        
    end
end

vtuFile = '../results/curvedBeamTspline';
msh_to_vtu (x, sigma, u, [noGPs noGPs], vtuFile);

%%
sigmaXX = sigma(1,:,:);
sigmaYY = sigma(2,:,:);
sigmaXY = sigma(3,:,:);
figure
scatter(xx(:,1),xx(:,2),160,sigmaXY(:),'full');
axis('equal');
colorbar






