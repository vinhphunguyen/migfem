%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Extended IGA for two dimensional multi-material elasticity problems.
%
% Plate with two materials.
% The enrichment function is the modified abs function proposed
% by Moes et al 2003.
%
% Vinh Phu Nguyen, April 2012
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

global p q controlPts W Q levelMax VOID CHI element weakEnrFunc


E1           = 1e2;  % Young modulus of matrix
nu1          = 0.2;  % Poisson ratio
E2           = 1e1;  % of inclusion
nu2          = 0.3;
stressState = 'PLANE_STRAIN';
L           = 1; % length of plate
levelMax    = 3;

vtuFile     = '../results/infinitePlateInclusionXIGA';

weakEnrFunc = 1; %1=Moes; 2=absolute of signed distance

plotEnrNodes=0; % do not plot interface and enriched nodes

% Elasticity matrix

Cm = elasticityMatrix(E1,nu1,stressState);
Ci = elasticityMatrix(E2,nu2,stressState);

tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%  CAD input: control points and knot vectors
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%lineInterfaceData
%lineInterfaceC2Data
lineInterfaceC0P2Data

noGPs          = p+1; % # of GPs, assume p>=q

% Boundary conditions

% 1. Bottom edge: fixed in y-direction
% 2. Top edge:    imposed y-displacement


% find boundary nodes for boundary conditions

bottomNodes = find(controlPts(:,2)==0)';
rightNodes  = find(controlPts(:,1)==L)';
leftNodes   = find(controlPts(:,1)==0)';
topNodes    = find(controlPts(:,2)==L)';

% essential boundary conditions

uNodes = bottomNodes(floor(noPtsX/2)+1);
vNodes = [bottomNodes topNodes];

uFixed     = zeros(size(uNodes));
vFixed     = [zeros(size(bottomNodes)) 0.1*ones(size(bottomNodes))];

udofs = 2*uNodes-1; % global indecies  of the fixed x disps
vdofs = 2*vNodes;   % global indecies  of the fixed y disps

% initialization

noDofs = noCtrPts*2 + size(inc_nodes,1)*1*2 ;

K = sparse(noDofs,noDofs);
f = zeros(noDofs,1);

% Due to the presence of additional dofs, the assembly is a little
% bit difficult than in FEM. We use fictitious nodes to handle these
% additional dofs. At an enriched node, we add one fantom node. These fictitious nodes
% are numbered from the total number of true nodes, ie, from numnode+1 ...

pos    = zeros(noCtrPts,1);
ninode = 0 ;

for i = 1 : noCtrPts
    if (enrich_node(i) == 3)
        pos(i) = (noCtrPts + ninode*1) + 1 ;
        ninode = ninode + 1 ;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Assembling system of equation
% Stiffness matrix and external force vector

coord=[-1 -1;1 -1;1 1;-1 1];

disp([num2str(toc),'  ASSEMBLING THE SYSTEM'])
%%
% Loop over elements (knot spans)

gps = [];

for e=1:noElems
    idu    = index(e,1);
    idv    = index(e,2);
    xiE    = elRangeU(idu,:); % [xi_i,xi_i+1]
    etaE   = elRangeV(idv,:); % [eta_j,eta_j+1]
    
    sctr   = element(e,:);         %  element scatter vector    
    pts    = controlPts(sctr,:);
    levelS = CHI(1,sctr);     
    
    W      = [];
    Q      = [];
    
    % -----------------------------------------------
    % Choose Gauss quadrature rules for elements
    
    voidId     = find(splitElems==e);
    if (isempty(voidId) == 0)     % split element
%         [aa] = hierarchicalGaussQuad(noGPs,CHI(elementV(e,:)),...
%             coord,node(elementV(e,:),:),splitElems(voidId,2),0);
        [W,Q] = discontQ4quad(7,chi(1,elementV(e,:)));
    else
        [W,Q] = quadrature(noGPs,'GAUSS',2);
    end
    
    % Determine the position in the global matrix K
    
    sctrB = assembly(e,enrich_node,pos);
    
    % loop over Gauss points
    
    for gp=1:size(W,1)
        pt      = Q(gp,:);
        wt      = W(gp);
        % compute coords in parameter space
        Xi      = parent2ParametricSpace(xiE, pt(1));
        Eta     = parent2ParametricSpace(etaE,pt(2));
        J2      = jacobianPaPaMapping(xiE,etaE);
        
        % compute derivative of basis functions w.r.t parameter coord
        
        [R dRdxi dRdeta] = NURBS2DBasisDers([Xi; Eta],p,q,uKnot,vKnot,weights');
        
        % compute the jacobian of physical and parameter domain mapping
        % then the derivative w.r.t spatial physical coordinates
        % Jacobian matrix
        
        levelset = dot(R,levelS);

        x        = R * pts; % global coord of GP
        gps      = [gps;x];
        %levelset = x(2)-0.5;
                
        if levelset >= 0
            C = Cm;
        else
            C = Ci;
        end
        
        % B matrix
        
        [B,J1] = BMatrixXIGA(Xi,Eta,e,enrich_node,R,dRdxi,dRdeta);
        
        % compute elementary stiffness matrix and
        % assemble it to the global matrix
        
        K(sctrB,sctrB) = K(sctrB,sctrB) + B' * C * B * J1 * J2 * wt;
    end
end


%%
disp([num2str(toc),'  APPLYING BOUNDARY  CONDITIONS'])

applyBC

% SOLVE SYSTEM

disp([num2str(toc),'  SOLVING THE SYSTEM'])

U = K\f;
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% POST-PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ux    = U(1:2:2*noCtrPts);
Uy    = U(2:2:2*noCtrPts);

%%

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
% exportfig(gcf,fileName,opts)

plot(gps(:,1),gps(:,2),'+');

% compute stresses and displacements at mesh vertices
disp([num2str(toc),'  COMPUTING STRESSES '])

%plot(gps(:,1),gps(:,2),'+');

plotStressXIGAMultiMats

% stressComp=4;
% figure
% clf
% plot_field(node,elementV,'Q4',Svm(:,1));
% hold on
% colorbar
% title('von Mises stress')
% axis off

fac = 5;

figure
clf
plot_field(node+fac*[dispX dispY],elementV,'Q4',dispX);
hold on
plot_mesh(node+fac*[dispX dispY],elementV,'Q4','b-');
colorbar
title('Displacement in x direction')
axis off

figure
clf
plot_field(node+fac*[dispX dispY],elementV,'Q4',dispY);
hold on
plot_mesh(node+fac*[dispX dispY],elementV,'Q4','b-');
colorbar
title('Displacement in y direction')
axis off

figure
clf
plot_field(node,elementV,'Q4',stress(:,:,3));
hold on
plot_mesh(node,elementV,'Q4','w-');
colorbar
title('Stress in y direction')

figure
clf
plot_field(node,elementV,'Q4',strain(:,:,3));
hold on
plot_mesh(node,elementV,'Q4','w-');
colorbar
title('Strain in x direction')

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,fileName,opts)








