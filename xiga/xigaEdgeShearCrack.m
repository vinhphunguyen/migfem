%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Extended Isogeometric analysis for two dimensional linear elastic
% fracture mechanics problems.
%
% Edge crack plate under shearing
%
% Vinh Phu Nguyen,
% Johns Hopkins University
% nvinhphu@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath ../fem_util/
addpath ../C_files/
addpath ../data/
addpath ../meshing/
addpath ../post-processing/
addpath ../fem-functions/
addpath ../analytical-solutions/
addpath ../nurbs-geopdes/inst/
addpath ../nurbs-util/

clc
clear all

global p q controlPts weights element xCrack xTips crack_node
global jDomainFac

E0          = 3e7;  % Young modulus
nu0         = 0.25;  % Poisson ratio
stressState ='PLANE_STRESS'; % either 'PLANE_STRAIN' or "PLANE_STRESS
sigma0      = 1;
jDomainFac  = 2;
vtuFile     = '../results/edgeCrackShear';

% COMPUTE ELASTICITY MATRIX

if ( strcmp(stressState,'PLANE_STRESS')  )      % Plane Strain case
    C=E0/(1-nu0^2)*[  1      nu0          0;
        nu0        1          0;
        0        0  (1-nu0)/2  ];
else                                            % Plane Strain case
    C=E0/(1+nu0)/(1-2*nu0)*[ 1-nu0      nu0        0;
        nu0    1-nu0        0;
        0        0  1/2-nu0 ];
end

tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%  CAD input: control points and knot vectors
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%edgeShearCrackC0Data 
%edgeShearCrackC2Data
edgeShearCrackCkData

noGPs         = p+1; % # of Gauss points along one direction
noGP1         = p+1;

% find boundary nodes for boundary conditions

fixedYNodes  =  find(controlPts(:,2)==0); %bottom edge
fixedXNodes  =  fixedYNodes;
forcedNodes  =  find(controlPts(:,2)==16); %top edge

% build connectivity ...

generateIGA2DMesh

% build boundary mesh for force vector computation

bndPoints      = controlPts(forcedNodes,:);
rightEdgeMesh  = zeros(noElemsU,p+1);

for i=1:noElemsU
    rightEdgeMesh(i,:) = forcedNodes(i:i+p);
end

% Each split node is enriched by ONE function, H(x)
% Each tip node is enriched by FOUR functions, B_i(x)
% then, the total dofs is :
% total dof = numnode*nsd + numsplitnode*1*nsd + numstipnode*4*nsd
% here, two dimension, nsd = 2

noDofs = noDofs + size(split_nodes,1)*1*2 + ...
    size(tip_nodes,1)*4*2;

K = sparse(noDofs,noDofs);  % global stiffness matrix
u = zeros(noDofs,1);        % displacement vector
f = zeros(noDofs,1);        % external force vector

% We use fictitious nodes/control points to handle these
% additional dofs. At a H(x) enriched node, we add one fantom node and
% at tip enriched node, four fantom nodes are added. These fictitious nodes
% are numbered from the total number of true nodes, ie, from numnode+1 ...

pos    = zeros(noCtrPts,1);
nsnode = 0 ;
ntnode = 0 ;

for i = 1 : noCtrPts
    if (enrich_node(i) == 1)
        pos(i) = (noCtrPts + nsnode*1 + ntnode*4) + 1 ;
        nsnode = nsnode + 1 ;
    elseif (enrich_node(i) == 2)
        pos(i) = (noCtrPts + nsnode*1 + ntnode*4) + 1 ;
        ntnode = ntnode + 1 ;
    end
end

jacob = zeros(2,2);

% essential boundary conditions

uFixed     = zeros(size(fixedXNodes))';
vFixed     = zeros(size(fixedYNodes))';

udofs      = 2*fixedXNodes-1;
vdofs      = 2*fixedYNodes;

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
    xiE    = elRangeU(idu,:); % [xi_i,xi_i+1]
    etaE   = elRangeV(idv,:); % [eta_j,eta_j+1]
    
    sctr   = element(e,:);          %  element scatter vector
    sctrV  = elementV(e,:);         %  element scatter vector Q4
    nn     = length(sctr);
    
    % -----------------------------------------------
    % Choose Gauss quadrature rules for elements
    
    if (ismember(e,split_elem))     % split element
        [W,Q] = discontQ4quad(7,levelSets(1,elementV(e,:),1));
    elseif (ismember(e,tip_elem))   % tip element
        [W,Q] = disTipQ4quad(7,levelSets(1,elementV(e,:)),node(sctrV,:),xTip);
    elseif ( any(intersect(tip_nodes,sctr)) ~= 0)% having tip enriched nodes
        [W,Q] = quadrature(10,'GAUSS',2);
    else
        [W,Q] = quadrature(noGPs,'GAUSS',2);
    end
    
    % Determine the position in the global matrix K
    
    sctrB = assembly(e,enrich_node,pos);
    
    % loop over Gauss points
    
    for gp=1:size(W,1)
        pt      = Q(gp,:);
      
        % compute coords in parameter space
        Xi      = parent2ParametricSpace(xiE,pt(1));
        Eta     = parent2ParametricSpace(etaE,pt(2));
        J2      = jacobianPaPaMapping(xiE,etaE);
    
        % compute derivative of basis functions w.r.t parameter coord

        [N dRdxi dRdeta] = NURBS2DBasisDers([Xi; Eta],p,q,uKnot,vKnot,weights');

        % B matrix
        
        [B,J1] = BMatrixXIGA(Xi,Eta,e,enrich_node,N,dRdxi,dRdeta);
        
        % Stiffness matrix
        
        K(sctrB,sctrB) = K(sctrB,sctrB) + B'*C*B*W(gp)*J1*J2;   
    end
end

% Computing external force

[W1,Q1] = quadrature(noGP1, 'GAUSS', 1 );

% Loop over elements along top edge 

for e=1:noElemsU
    xiE   = elRangeU(e,:); % [xi_i,xi_i+1]
    conn  = elConnU(e,:);
    pts   = bndPoints(conn,:);
    
    sctrx = 2*rightEdgeMesh(e,:)-1;
    
    % loop over Gauss points
    for gp=1:size(W1,1)
        xi      = Q1(gp,:);
        wt      = W1(gp);
        Xi      = 0.5 * ( ( xiE(2) - xiE(1) ) * xi + xiE(2) + xiE(1));
        J2      = 0.5 * ( xiE(2) - xiE(1) );
        
        % compute derivative of basis functions w.r.t parameter coord
        
        [N dNdxi] = NURBS1DBasisDers(Xi,p,uKnot,weights); 
        
        % compute the jacobian of physical and parameter domain mapping
        % then the derivative w.r.t spatial physical coordinates
        
        jacob1   = dNdxi * pts;
        J1       = norm (jacob1);
              
        f(sctrx) = f(sctrx) + N' * sigma0 * J1 * J2 * wt;
    end
end

disp([num2str(toc),'  APPLYING BOUNDARY  CONDITIONS'])
bcwt=mean(diag(K)); % a measure of the average  size of an element in K
% used to keep the  conditioning of the K matrix

applyBC


% SOLVE SYSTEM
disp([num2str(toc),'  SOLVING THE SYSTEM'])
U=K\f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% POST-PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ux = U(1:2:noDofs);
Uy = U(2:2:noDofs);

% compute stresses and displacements at mesh vertices
disp([num2str(toc),'  COMPUTING STRESSES '])

plotStressXIGA

% plot stress

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
title('Displacement in y direction')
axis off

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,fileName,opts)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp([num2str(toc),'      SIFs COMPUTATION'])

computeInteractionIntegrals

% Compute SIFs from I integral

if ( strcmp(stressState,'PLANE_STRAIN') )
    Eb = E0/(1-nu0^2);
else
    Eb = E0;
end

Knum   = 0.5*Eb*I

% Exact SIFs from Yau Wang Corten 1980

Kexact = [34 4.55]

% Relative error

Knum(1)/Kexact(1)
Knum(2)/Kexact(2)

%%%%%%%%%%%%

% crLength = [0.4; 0.45];
% KI       = [2.3757;2.8975];
% KIExact  = [2.3585;2.8772];
% hold on
% cr1 = plot(crLength,KI,'-');
% cr2 = plot(crLength,KIExact,'*');


