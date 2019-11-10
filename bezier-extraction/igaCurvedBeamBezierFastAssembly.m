%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Isogeometric analysis for two dimensional elasticity problems.
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

noGPs         = 3; % # of Gauss points along one direction


E0          = 3e7;  % Young modulus
nu0         = 0.25;  % Poisson ratio
stressState ='PLANE_STRESS'; % either 'PLANE_STRAIN' or "PLANE_STRESS
pressure    = 3e4; % inner pressure

% COMPUTE ELASTICITY MATRIX

De = elasticityMatrix(E0,nu0,stressState);

tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%  CAD input: control points and knot vectors
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%annularDataGeopdesBezier
annularDataP2Q1Geopdes % not the same file in ../data/

% find boundary nodes for boundary conditions

fixedXNodes  =  find(controlPts(:,1)==0);
fixedYNodes  =  find(controlPts(:,2)==0);
forcedNodes  =  1:noPtsX:noCtrPts;

% build connectivity ...

generateIGA2DMesh

% build boundary mesh for force vector computation

bndPoints      = controlPts(forcedNodes,:);
rightEdgeMesh  = zeros(noElemsV,q+1);

for i=1:noElemsV
    rightEdgeMesh(i,:) = forcedNodes(i:i+q);
end

% initialization

u = zeros(noDofs,1);        % displacement vector
f = zeros(noDofs,1);        % external force vector

%% one the assembly
nElNod = size(element,2);
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

udofs      = fixedXNodes;
vdofs      = fixedYNodes+noCtrPts;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gauss quadrature rule
[W,Q]=quadrature(  noGPs, 'GAUSS', 2 ); % noGPs x noGPs point quadrature

%% Pre-compute Bernstein basis and derivatives for ONE Bezier element

noBasis = (p+1)*(q+1);
noGpEle = noGPs*noGPs;

shapes  = zeros(noGpEle,noBasis);
derivs  = zeros(noGpEle,noBasis,2);

for gp=1:size(W,1)
    [shapes(gp,:) derivs(gp,:,:)] = getShapeGrad2Bernstein2D(p,q,Q(gp,1),Q(gp,2));
end

%% Assembling system of equation
% Stiffness matrix and external force vector

disp([num2str(toc),'  ASSEMBLING THE SYSTEM'])

% Loop over elements 

for e=1:noElems
    sctr   = element(e,:);         %  element scatter vector
    sctrB  = [sctr sctr+noCtrPts]; %  vector that scatters a B matrix
    nn     = length(sctr);
    
    B      = zeros(3,2*nn);
    Ke     = Ke0;
    
    Ce     = C(:,:,e);             % element Bezier extraction operator
    we     = diag(weights(sctr,:));% element weights
    pts    = controlPts(sctr,:);   % element nodes 
    Wb     = Ce'*weights(sctr,:);  % element Bezier weights
    
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
    
    jSprRw(jSprLK) = mRwGrd(:)
    jSprCl(jSprLK) = mClGrd(:)
    vSprGK(jSprLK) = Ke(:);
    
    jSprLK         = jSprLK + nElmLK; % move to the next position 
    
end


% Here comes the total stiffness matrix in one shot!!!

K = sparse(jSprRw,jSprCl,vSprGK,noDofs,noDofs);

%% Computing external force

[W1,Q1] = quadrature(2, 'GAUSS', 1 );

% Loop over elements along left edge = noElemsV
% which is a NURBS curve of order q

for e=1:noElemsV    
    sctrx = rightEdgeMesh(e,:);
    sctry = sctrx + noCtrPts;
    
    CetE  = Cet(:,:,e);              % Bezier extration operator
    we    = diag(weights(sctrx,:));  % element weights
    Wb    = CetE'*weights(sctrx,:);  % element Bezier weights
    pts   = controlPts(sctrx,:);
    
    % loop over Gauss points
    for gp=1:size(W1,1)
        xi      = Q1(gp,:);
        wt      = W1(gp);
     
        [Be,dBdxi]= bernsteinBasis1D(q,xi);        
        wb        = dot(Be,Wb);        
        R         = we*CetE*Be/wb;
        dwbdxi    = dot(dBdxi,Wb);
        dRdxi     = we*CetE*(dBdxi/wb-dwbdxi*Be/(wb*wb));
        
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

%%
disp([num2str(toc),'  APPLYING BOUNDARY  CONDITIONS'])
bcwt=mean(diag(K)); % a measure of the average  size of an element in K
% used to keep the  conditioning of the K matrix

applyBC


% SOLVE SYSTEM
disp([num2str(toc),'  SOLVING THE SYSTEM'])
U=K\f;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% POST-PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% There are three techniques:
% 1. Using Gauss point information and a triangulation of GPs
% 2. Using the so-called visualization Q4 mesh (requires knots information)
% 3. Using Gauss point information and a mesh of Q4 elements for each NURBS
%    elements. There are no connection between Q4 elements of neighbouring
%    NURBS elements. Therefore, there are empty spaces in the
%    visualization. However implementation is easy and works for
%    NURBS/Tsplines using Bezier extraction.

Ux    = U(1:noCtrPts);
Uy    = U(noCtrPts+1:noDofs);

vtuFile     = '../results/curvedBeamNURBSFA';


buildVisualizationMesh;
stress = zeros(noElems,size(elementV,2),3);
disp   = zeros(noElems,size(elementV,2),2);

for e=1:noElems
    idu    = index(e,1);
    idv    = index(e,2);
    xiE    = elRangeU(idu,:); % [xi_i,xi_i+1]
    etaE   = elRangeV(idv,:); % [eta_j,eta_j+1]
    
    sctr   = element(e,:);         %  element scatter vector
    sctrB  = [sctr sctr+noCtrPts]; %  vector that scatters a B matrix
    nn     = length(sctr);
    
    B      = zeros(3,2*nn);
    pts    = controlPts(sctr,:);
    
    uspan = FindSpan(noPtsX-1,p,xiE(1), uKnot);
    vspan = FindSpan(noPtsY-1,q,etaE(1),vKnot);
    
    % loop over Gauss points
    
    gp = 1;
    for iv=1:2
        if (iv==2)
            xiE = sort(xiE,'descend');
        end
        for iu=1:2            
            Xi  = xiE(iu);
            Eta = etaE(iv);
            [N dRdxi dRdeta] = NURBS2DBasisDersSpecial([Xi; Eta],...
                p,q,uKnot,vKnot,weights',[uspan;vspan]);
            
            % compute the jacobian of physical and parameter domain mapping
            % then the derivative w.r.t spatial physical coordinates
            
            jacob      = pts' * [dRdxi' dRdeta'];
            
            if (det(jacob) <= 1e-8)
                [N dRdxi dRdeta] = NURBS2DBasisDers([Xi-0.01; ...
                    Eta],p,q,uKnot,vKnot,weights');
                jacob      = pts' * [dRdxi' dRdeta'];
            end
            
            % Jacobian inverse and spatial derivatives
            
            invJacob   = inv(jacob);
            dRdx       = [dRdxi' dRdeta'] * invJacob;
            
            % B matrix
            
            B(1,1:nn)       = dRdx(:,1)';
            B(2,nn+1:2*nn)  = dRdx(:,2)';
            B(3,1:nn)       = dRdx(:,2)';
            B(3,nn+1:2*nn)  = dRdx(:,1)';
            
            strain          = B*U(sctrB);
            stress(e,gp,:)  = De*strain;
            disp(e,gp,:)    = N*[Ux(sctr) Uy(sctr)];
            
            gp = gp +1;
        end
    end
end

stressComp=3;
figure
clf
plot_field(node,elementV,'Q4',stress(:,:,stressComp));
hold on
colorbar
title('Stress in x direction')
axis off
%plot_mesh(node,elementV,'Q4','g.-');

figure
clf
plot_field(node,elementV,'Q4',disp(:,:,2));
hold on
colorbar
title('Displacement in x direction')
axis off

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,fileName,opts)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% export to VTK format to plot in Mayavi or Paraview

sigmaXX = zeros(size(node,1),1);
sigmaYY = zeros(size(node,1),1);
sigmaXY = zeros(size(node,1),1);

dispX = zeros(size(node,1),1);
dispY = zeros(size(node,1),1);

for e=1:size(elementV,1)
    connect = elementV(e,:);
    for in=1:4
        nid = connect(in);
        sigmaXX(nid) = stress(e,in,1);
        sigmaYY(nid) = stress(e,in,2);
        sigmaXY(nid) = stress(e,in,3);
        
        dispX(nid) = disp(e,in,1);
        dispY(nid) = disp(e,in,2);
    end
end

VTKPostProcess(node,elementV,2,'Quad4',vtuFile,...
             [sigmaXX sigmaYY sigmaXY],[dispX dispY]);



%% Visualizaition techqique 3

dRdxi=[];
x     = zeros(3,noGpEle,noElems);  % global coords of Gauss points
u     = zeros(3,noGpEle,noElems);  % displacements of Gauss points
sigma = zeros(3,noGpEle,noElems);  % stresses      of Gauss points

for e=1:noElems
    sctr   = element(e,:);         %  element scatter vector
    sctrB  = [sctr sctr+noCtrPts]; %  vector that scatters a B matrix
    nn     = length(sctr);
    
    B      = zeros(3,2*nn);
    Ke     = Ke0;
    
    Ce     = C(:,:,e);                 % element Bezier extraction operator
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
        u(1:2,gp,e)     = R' * [Ux(sctr) Uy(sctr)];
        sigma(:,gp,e)   = De * B * elemDisp;
                
    end
end

msh_to_vtu (x, sigma, u, [noGPs noGPs], vtuFile);








