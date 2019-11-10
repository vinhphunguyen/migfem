%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Isogeometric analysis for two dimensional elasticity problems.
% Tspline geometry.
% Using the so-called Bezier extraction operator.
% Using the triple sparse matrix format to speed up the assembly process.
% This works for cases in which elements have the same number of nodes.
%
% A square specimen in tension. Serve to verify the implementation.
% Both structured and unstructured Tspline input files are illustrated.
%
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

E0          = 250;           % Young modulus
nu0         = 0.3;           % Poisson ratio
stressState ='PLANE_STRAIN'; % either 'PLANE_STRAIN' or "PLANE_STRESS
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

%filename = 'square_structured.iga';
filename = 'square_unstructured.iga';
tspline  = read_bezier_extraction (filename);
convert2DTsplines;


%% initialization

K = sparse(noDofs,noDofs);  % stiffness matrix
u = zeros(noDofs,1);        % displacement vector
f = zeros(noDofs,1);        % external force vector


%% essential boundary conditions

EPS = 1e-8;
fixedNodes  = find(abs(controlPts(:,1))<EPS)';
bottomNodes = find(abs(controlPts(:,2))<EPS)';
forcedNodes = find(abs(controlPts(:,1)-1)<EPS)';

uFixed     = zeros(size(fixedNodes));
vFixed     = zeros(1);

udofs=fixedNodes;          % global indecies  of the fixed x displacements
vdofs=intersect(fixedNodes,bottomNodes)+noCtrPts;


uFixed = [uFixed ubar*ones(1,length(forcedNodes))];
udofs = [udofs forcedNodes];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Assembling system of equation
% Stiffness matrix and external force vector

disp([num2str(toc),'  ASSEMBLING THE SYSTEM'])

% Loop over elements

for e=1:noElems
    sctr   = element{e};           %  element scatter vector
    sctrB  = [sctr sctr+noCtrPts]; %  vector that scatters a B matrix
    nn     = length(sctr);    
    B      = zeros(3,2*nn);
    
    degree = degrees(e,:);         % Bernstein basis orders in two dirs.
    Ce     = C{e};                 % element Bezier extraction operator
    we     = diag(weights(sctr));  % element weights
    pts    = controlPts(sctr,:);   % element nodes
    Wb     = Ce'*weights(sctr);    % element Bezier weights
    
    % retrieve correctly the Bernstein basis
    % of the current element
    bid    = find(ismember(uDegree,degree),1);
    shapes = bernstein(bid).basis;
    derivs = bernstein(bid).ders;    
    W      = quad(bid).weights;
        
    % loop over Gauss points
    
    for gp=1:size(W,1)        
        wt      = W(gp);
        
        %% Bernstein basis and derivatives at GP gp
        Be      = shapes(gp,:)';
        dBedxi  = reshape(derivs(gp,:,:),prod(degree+1),2);

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
        
        K(sctrB,sctrB) = K(sctrB,sctrB) + B' * De * B * detJ * wt;
    end
    dRdxi=[];
end

%% Computing external force

%f(forcedNodes+noCtrPts)=f0;

%%
disp([num2str(toc),'  APPLYING BOUNDARY  CONDITIONS'])
bcwt=mean(diag(K)); % a measure of the average  size of an element in K
% used to keep the  conditioning of the K matrix

applyBC


% SOLVE SYSTEM
disp([num2str(toc),'  SOLVING THE SYSTEM'])
U=K\f;


%%
%%%
%%% POST-PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ux    = U(1:noCtrPts);
Uy    = U(noCtrPts+1:noDofs);

x     = zeros(3,noGpEle,noElems);  % global coords of Gauss points
u     = zeros(3,noGpEle,noElems);  % displacements of Gauss points
sigma = zeros(3,noGpEle,noElems);  % stresses at Gauss points

ind = 1;

for e=1:noElems
    sctr   = element{e};         %  element scatter vector
    sctrB  = [sctr sctr+noCtrPts]; %  vector that scatters a B matrix
    nn     = length(sctr);
    
    B      = zeros(3,2*nn);
    
    degree = degrees(e,:);
    Ce     = C{e};                 % element Bezier extraction operator
    we     = diag(weights(sctr));  % element weights
    pts    = controlPts(sctr,:);   % element nodes
    Wb     = Ce'*weights(sctr);    % element Bezier weights
    
    elemDisp = U(sctrB);
    
    % retrieve correctly the Bernstein basis
    % of the current element
    bid    = find(ismember(uDegree,degree),1);
    shapes = bernstein(bid).basis;
    derivs = bernstein(bid).ders;
    
    W      = quad(bid).weights;
    Q      = quad(bid).points;
    
    % loop over Gauss points
    
    for gp=1:size(W,1)
        pt      = Q(gp,:);
        wt      = W(gp);
        
        %% Bernstein basis and derivatives at GP gp
        Be      = shapes(gp,:)';
        dBedxi  = reshape(derivs(gp,:,:),prod(degree+1),2);
        
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
        u(1:2,gp,e)     = R' * [Ux(sctr) Uy(sctr)];
        sigma(:,gp,e)   = De * B * elemDisp;
        
        ind = ind + 1;
    end
    dRdxi=[];
end

msh_to_vtu (x, sigma,u,[max(max(uDegree))+1 max(max(uDegree))+1], '../results/squareTspline');





