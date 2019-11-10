%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Isogeometric analysis for two dimensional elasticity problems.
% Using the so-called Bezier extraction operator.
%
% Curved beam subject to end shear. 
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

fixedXNodes1  =  find(controlPts(:,1)==0);
fixedYNodes   =  fixedXNodes1;

fixedXNodes2  =  find(controlPts(:,2)==0);

fixedXNodes   = [fixedXNodes1; fixedXNodes2];

% build connectivity ...

generateIGA2DMesh


% initialization

K = sparse(noDofs,noDofs);  % global stiffness matrix
u = zeros(noDofs,1);        % displacement vector
f = zeros(noDofs,1);        % external force vector

% essential boundary conditions

uFixed     = [zeros(size(fixedXNodes1))' ones(size(fixedXNodes2))'];
vFixed     = zeros(size(fixedYNodes))';

udofs      = 2*fixedXNodes-1;
vdofs      = 2*fixedYNodes;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gauss quadrature rule

noGPs         = p+1; % # of Gauss points along one direction

[W,Q]=quadrature( noGPs, 'GAUSS', 2 ); % noGPs x noGPs point quadrature

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
    sctr   = element(e,:);         %  element scatter vector
    nn     = length(sctr);
    
    sctrB(1,1:2:2*nn) = 2*sctr-1;
    sctrB(1,2:2:2*nn) = 2*sctr;
    
    B      = zeros(3,2*nn);
    
    Ce     = C(:,:,e);             % element Bezier extraction operator
    we     = diag(weights(sctr));  % element weights
    pts    = controlPts(sctr,:);   % element nodes
    Wb     = Ce'*weights(sctr);    % element Bezier weights
    
    % loop over Gauss points
    
    for gp=1:size(W,1)
        pt      = Q(gp,:);
        wt      = W(gp);
        
        %% Bernstein basis and derivatives at GP gp
        Be      = shapes(gp,:)';
        dBedxi  = reshape(derivs(gp,:,:),noBasis,2);
        
        %% Bezier weight functions (denomenator of NURBS)
        wb        = dot(Be,Wb);            % Be(I)*Wb(I)
        dwbdxi(1) = dot(dBedxi(:,1),Wb);   % Be(I)_{,xi} * Wb(I)
        dwbdxi(2) = dot(dBedxi(:,2),Wb);   % Be(I)_{,et} * Wb(I)
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
        
        K(sctrB,sctrB) = K(sctrB,sctrB) + B' * De * B * detJ * wt;
    end
end

disp([num2str(toc),'  APPLYING BOUNDARY  CONDITIONS'])

[K,f]=applyDirichletBCs(K,f,udofs,vdofs,uFixed,vFixed);

% SOLVE SYSTEM
disp([num2str(toc),'  SOLVING THE SYSTEM'])
U=K\f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% POST-PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ux    = U(1:2:noDofs);
Uy    = U(2:2:noDofs);

vtuFile     = '../results/curvedBeamShearNURBS';


buildVisualizationMesh;
stress = zeros(noElems,size(elementV,2),3);
disp   = zeros(noElems,size(elementV,2),2);

for e=1:noElems
    idu    = index(e,1);
    idv    = index(e,2);
    xiE    = elRangeU(idu,:); % [xi_i,xi_i+1]
    etaE   = elRangeV(idv,:); % [eta_j,eta_j+1]
    
    sctr   = element(e,:);         %  element scatter vector
    nn     = length(sctr);
    sctrB(1,1:2:2*nn) = 2*sctr-1;
    sctrB(1,2:2:2*nn) = 2*sctr;
    
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
            B(1,1:2:2*nn)  = dRdx(:,1)';
            B(2,2:2:2*nn)  = dRdx(:,2)';
            B(3,1:2:2*nn)  = dRdx(:,2)';
            B(3,2:2:2*nn)  = dRdx(:,1)';
            
            strain          = B*U(sctrB);
            stress(e,gp,:)  = De*strain;
            disp(e,gp,:)    = N*[Ux(sctr) Uy(sctr)];
            
            gp = gp +1;
        end
    end
end

stressComp=1;
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











