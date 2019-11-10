%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Isogeometric analysis for two dimensional elasticity problems.
% Tspline geometry.
% Using the so-called Bezier extraction operator.
% Using the triple sparse matrix format to speed up the assembly process.
% 
% Procedure:
%
% 1. Using Rhino3d and Tspline plugin to make the geometry. Save it as a
% *.iga file.
%
% 2. Run IGA on that file.
%
% A plate with a centered circular hole in tension.
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


E0          = 1e5;           % Young modulus
nu0         = 0.3;           % Poisson ratio
stressState ='PLANE_STRESS'; % either 'PLANE_STRAIN' or "PLANE_STRESS
f0          = 10.;
ubar=.4;

% COMPUTE ELASTICITY MATRIX

De = elasticityMatrix(E0,nu0,stressState);

tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%  CAD input: *.iga file from Rhino3d plugin
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filename = 'plate_with_hole.iga';
tspline  = read_bezier_extraction (filename);
convert2DTsplines;


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

EPS = 1e-8;
fixedXNodes  = find(abs(controlPts(:,1))<EPS)';
fixedYNodes  = find(abs(controlPts(:,2))<EPS)';
forcedNodes  = find(abs(controlPts(:,2)-4)<EPS)'; % top edge

uFixed     = zeros(size(fixedXNodes));
vFixed     = zeros(size(fixedYNodes));

udofs=fixedXNodes;          % global indecies  of the fixed x displacements
vdofs=fixedYNodes+noCtrPts; % global indecies  of the fixed y displacements

% non-zero Dirichlet BCs

vFixed = [vFixed ubar*ones(1,length(forcedNodes))];
vdofs  = [vdofs forcedNodes+noCtrPts];

%% check BCs

figure
hold on
plot(controlPts(:,1),controlPts(:,2),'r*');
axis('equal');
plot(controlPts(fixedXNodes,1),controlPts(fixedXNodes,2),'cs','MarkerSize',12);
plot(controlPts(fixedYNodes,1),controlPts(fixedYNodes,2),'cs','MarkerSize',12);
plot(controlPts(forcedNodes,1),controlPts(forcedNodes,2),'bo','MarkerSize',12);

%%

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
    sctr   = element{e};         %  element scatter vector
    sctrB  = [sctr sctr+noCtrPts]; %  vector that scatters a B matrix
    nn     = length(sctr);
    
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
        B(1,1:nn)       = dRdx(:,1)';
        B(2,nn+1:2*nn)  = dRdx(:,2)';
        B(3,1:nn)       = dRdx(:,2)';
        B(3,nn+1:2*nn)  = dRdx(:,1)';
        
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

%f(forcedNodes+noCtrPts)=f0;

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

disp([num2str(toc),'  POST-PROCESSING'])

Ux    = U(1:noCtrPts);
Uy    = U(noCtrPts+1:noDofs);

xx    = zeros(noGpEle*noElems,2);  % global coords of Gauss points
x     = zeros(3,noGpEle,noElems);  % global coords of Gauss points
u     = zeros(3,noGpEle,noElems);  % global coords of Gauss points
sigma = zeros(3,noGpEle,noElems);  % stresses at Gauss points

id=1;
for e=1:noElems
    sctr   = element{e};         %  element scatter vector
    sctrB  = [sctr sctr+noCtrPts]; %  vector that scatters a B matrix
    nn     = length(sctr);
    
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
        B(1,1:nn)       = dRdx(:,1)';
        B(2,nn+1:2*nn)  = dRdx(:,2)';
        B(3,1:nn)       = dRdx(:,2)';
        B(3,nn+1:2*nn)  = dRdx(:,1)';
        
        x(1:2,gp,e)     = R' * pts;
        xx(id,:)        = R' * pts;
        u(1:2,gp,e)     = R' * [Ux(sctr) Uy(sctr)];
        sigma(:,gp,e)   = De * B * elemDisp;
       id = id + 1;         
    end
end

tri = delaunay(xx);

r_in=1;
r_out=4;
numy=10;
numx=10;
count_tri=size(tri,1);
i=1;count=1;aa=1;
while count <= count_tri
    cord = xx (tri(i,:),:);
    if max(max((cord)))<r_in+(r_out-r_in)/numy   %
        for j=1:size(cord,1)
            if j == size(cord,1)
                x1=cord (j,1); y1=cord(j,2);
                x2=cord(1,1) ; y2=cord(1,2);
            else
                x1=cord (j,1)  ; y1=cord(j,2);
                x2=cord(j+1,1) ; y2=cord(j+1,2);
            end
            ll(j) = sqrt((x2-x1)^2+(y1-y2)^2);
        end
        
        %if max(ll) >= 1.
        if max(ll) >= r_out*pi/numx/(noGPs-1)
            tri(i,:)=[];
            i=i-1;
        end
    end
    count = count+1;
    i=i+1;
end

sigmaX = sigma(1,:,:);
sigmaY = sigma(2,:,:);
ux     = u(1,:,:);
uy     = u(2,:,:);

figure
hold on
plot_field(xx,tri,'T3',sigmaY(:));
axis('equal');
xlabel('X');
ylabel('Y');
title('Sigma XX');
set(gcf,'color','white');
colorbar('vert');
%triplot(tri)

figure
hold on
plot_field(xx,tri,'T3',ux(:));
axis('equal');
xlabel('X');
ylabel('Y');
title('ux');
set(gcf,'color','white');
colorbar('vert');

figure
hold on
plot_field(xx,tri,'T3',uy(:));
axis('equal');
xlabel('X');
ylabel('Y');
title('uy');
set(gcf,'color','white');
colorbar('vert');

% discrete stress contour
% if there are a lots of GPs then it works well
figure
scatter(xx(:,1),xx(:,2),80,sigmaY(:),'full');
axis('equal');
colorbar

%
vtuFile     = '../results/plateHoleTspline';
msh_to_vtu (x, sigma, u, [noGPs noGPs], vtuFile);
