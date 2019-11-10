%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Extended Isogeometric analysis for two dimensional linear elastic
% fracture mechanics problems.
%
% Center crack infinite plate in tension
% Dirichlet BCs = exact displacements imposed with the
% Least Square method (Luycker et al. IJNME 2011).
%
% Vinh Phu Nguyen,
% Delft University of Technology/Johns Hopkins University
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
global uKnot vKnot noDofs levelSets split_nodes spit_elem

E0          = 1e6;  % Young modulus
nu0         = 0.3;  % Poisson ratio
stressState ='PLANE_STRAIN'; % either 'PLANE_STRAIN' or "PLANE_STRESS
sigmato     = 1e4;

vtuFile      = 'infinieCrackModeILS';
vtuCrackFile = 'infinieCrackModeILS-cracked';

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

edgeInfiniteCenterCrackC0Data
%edgeInfiniteCenterCrackC2Data
%edgeInfiniteCenterCrackC4Data

noGPs         = p+1; % # of Gauss points along one direction
noGP1         = p+1;

% find boundary nodes for boundary conditions

bottomNodes =  find(node(:,2)==0);
rightNodes  =  find(node(:,1)==D);
topNodes    =  find(node(:,2)==D);
leftNodes   =  find(node(:,1)==0);

bottomNodesIGA =  find(controlPts(:,2)==0);
rightNodesIGA  =  find(controlPts(:,1)==D);
topNodesIGA    =  find(controlPts(:,2)==D);
leftNodesIGA   =  find(controlPts(:,1)==0);

topNodesIGA = sort(topNodesIGA,'descend');

% Dirichlet control points/nodes
% should not use function "unique" because it will
% automatically sort the vector!!!

dispNodes   = [bottomNodesIGA ; rightNodesIGA(2:end);...
               topNodesIGA(2:end)];
noDispNodes = length(dispNodes);

% build boundary mesh

bottomPoints    = controlPts(bottomNodesIGA,:);
rightPoints     = controlPts(rightNodesIGA,:);
topPoints       = controlPts(topNodesIGA,:);
leftPoints      = controlPts(leftNodesIGA,:);

% linear two node line elements

bottomEdgeMesh  = zeros(length(bottomNodesIGA)-1,2);
rightEdgeMesh   = zeros(length(rightNodesIGA)-1,2);
topEdgeMesh     = zeros(length(topNodesIGA)-1,2);
leftEdgeMesh    = zeros(length(leftNodesIGA)-1,2);

bottomEdgeMeshIGA  = zeros(noElemsU,p+1);
rightEdgeMeshIGA   = zeros(noElemsV,q+1);
topEdgeMeshIGA     = zeros(noElemsU,p+1);
leftEdgeMeshIGA    = zeros(noElemsU,q+1);

for i=1:noElemsU
    bottomEdgeMeshIGA(i,:) = bottomNodesIGA(i:i+p);
    topEdgeMeshIGA(i,:)    = topNodesIGA(i:i+p);
end

for i=1:noElemsV
    rightEdgeMeshIGA(i,:) = rightNodesIGA(i:i+p);
    leftEdgeMeshIGA(i,:)  = leftNodesIGA(i:i+p);
end

% for i=1:noElemsU
%     bottomEdgeMesh(i,:) = bottomNodes(i:i+1);
%     topEdgeMesh(i,:)    = topNodes(i:i+1);
% end
%
% for i=1:noElemsV
%     rightEdgeMesh(i,:) = rightNodes(i:i+1);
%     leftEdgeMesh(i,:)  = leftNodes(i:i+1);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Least square method for determining the control
% points displacements

disp([num2str(toc),'  LEAST SQUARE for DIRICHLET BCs'])

noBndElems = noElemsU*2 + noElemsV;
bndElement = zeros(noBndElems,p+1);

for ie=1:noElemsU
    bndElement(ie,:) = [ie:ie+p];
end

start = bottomNodesIGA(end);

for ie=1:noElemsV
    bndElement(ie+noElemsU,:) = [start:start+p];
    start = start + 1;
end

start = max(max(bndElement));

for ie=1:noElemsU
    bndElement(ie+noElemsU+noElemsV,:) = [start:start+p];
    start = start + 1;
end

A  = zeros(noDispNodes,noDispNodes);
bx = zeros(noDispNodes,1);
by = zeros(noDispNodes,1);

% bottom edge

noxC   = 3;

for ie=1:noElemsU
    sctr   = bottomEdgeMeshIGA(ie,:);
    pts    = controlPts(sctr,:);
    sctrA  = bndElement(ie,:);
    xiE    = elRangeU(ie,:);
    xiArr  = linspace(xiE(1),xiE(2),noxC);
   
    for ic=1:noxC                
        xi = xiArr(ic);
        [N dNdxi] = NURBS1DBasisDers(xi,p,uKnot,weights);
        
        A(sctrA,sctrA) = A(sctrA,sctrA) + N'*N;
        
        x        = N    *pts; 
        
        % exact displacements
        
        [ux,uy] = exactDispModeI(x,E0,nu0,stressState,sigmato,xTip,seg,...
            cracklength);
        
        bx(sctrA) = bx(sctrA) + ux*N';
        by(sctrA) = by(sctrA) + uy*N';
    end
end

% right edge

for ie=1:noElemsV
    sctr   = rightEdgeMeshIGA(ie,:);
    pts    = controlPts(sctr,:);
    sctrA  = bndElement(ie+noElemsU,:);
    xiE    = elRangeV(ie,:);
    xiArr  = linspace(xiE(1),xiE(2),noxC);
    
    for ic=1:noxC                      
        xi             = xiArr(ic);        
        [N dNdxi]      = NURBS1DBasisDers(xi,q,vKnot,weights);        
        A(sctrA,sctrA) = A(sctrA,sctrA) + N'*N;        
        x              = N*pts; 
        
        % exact displacements
        
        [ux,uy] = exactDispModeI(x,E0,nu0,stressState,...
                                 sigmato,xTip,seg,cracklength);
        
        bx(sctrA) = bx(sctrA) + ux*N';
        by(sctrA) = by(sctrA) + uy*N';
    end
end

% top edge

for ie=1:noElemsU
    sctr   = topEdgeMeshIGA(ie,:);
    pts    = controlPts(sctr,:);    
    sctrA  = bndElement(ie+noElemsU+noElemsV,:);
    xiE    = elRangeU(ie,:);
    xiArr  = linspace(xiE(1),xiE(2),noxC);
        
    for ic=1:noxC                        
        xi        = xiArr(ic);        
        [N dNdxi] = NURBS1DBasisDers(xi,p,uKnot,weights);        
        A(sctrA,sctrA) = A(sctrA,sctrA) + N'*N;        
        x         = N*pts; 
        
        % exact displacements
        
        [ux,uy] = exactDispModeI(x,E0,nu0,stressState,...
                                 sigmato,xTip,seg,cracklength);
        
        bx(sctrA) = bx(sctrA) + ux*N';
        by(sctrA) = by(sctrA) + uy*N';
    end
end

% solve the system Aq_x=bx and Aq_y=by

[LL UU] = lu(A);
qxTemp  = LL\bx;
qyTemp  = LL\by;
qx      = UU\qxTemp;
qy      = UU\qyTemp;

%%%%

uFixed     = qx';
vFixed     = qy';

udofs      = 2*dispNodes-1;
vdofs      = 2*dispNodes;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    nn     = length(sctr);
    
    % -----------------------------------------------
    % Choose Gauss quadrature rules for elements
    
    if     (ismember(e,split_elem))              % split element
        [W,Q] = quadrature(12,'GAUSS',2);
    elseif (ismember(e,tip_elem))                % tip element
        [W,Q] = quadrature(12,'GAUSS',2);
    elseif ( any(intersect(tip_nodes,sctr)) ~= 0)% having tip enriched nodes
        [W,Q] = quadrature(12,'GAUSS',2);
    else
        [W,Q] = quadrature(noGPs,'GAUSS',2);
    end
    
    % Determine the position in the global matrix K
    
    sctrB = assembly(e,element,enrich_node,pos);
    
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

disp([num2str(toc),'  APPLYING BOUNDARY  CONDITIONS'])

[W1,Q1] = quadrature(5, 'GAUSS', 1 );

% Loop over elements along left edge = noElemsV
% left edge: natural BCs

for e=1:noElemsV
    xiE   = elRangeV(e,:); % [xi_i,xi_i+1]
    conn  = elConnV(e,:);
    sctrx = 2*leftEdgeMeshIGA(e,:)-1;
    sctry = 2*leftEdgeMeshIGA(e,:);
    pts   = leftPoints(conn,:);
    
    % loop over Gauss points
    for gp=1:size(W1,1)
        xi      = Q1(gp,:);
        wt      = W1(gp);
        Xi      = 0.5 * ( ( xiE(2) - xiE(1) ) * xi + xiE(2) + xiE(1));
        J2      = 0.5 * ( xiE(2) - xiE(1) );
        
        [N dNdxi] = NURBS1DBasisDers(Xi,p,uKnot,weights);
        
        % compute the jacobian of physical and parameter domain mapping
        % then the derivative w.r.t spatial physical coordinates
        
        x        = N    *pts; % global coord of GP
        jacob1   = dNdxi*pts;
        J1       = norm (jacob1);
        
        % compute exact stresses
        
        [sigmaxx,sigmayy,sigmaxy] = exactStressModeI(x,E0,nu0,...
            sigmato,xTip,seg,cracklength);
        
        % then the tractions, note that normal vector is [-1,0]
        
        tx = -sigmaxx;
        ty = -sigmaxy;
        
        f(sctrx) = f(sctrx) + N' * tx * J1 * J2 * wt;
        f(sctry) = f(sctry) + N' * ty * J1 * J2 * wt;
    end
end

% SOLVE SYSTEM

disp([num2str(toc),'  SOLVING THE SYSTEM'])

applyBC

U = K\f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% POST-PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ux = U(1:2:noDofs);
Uy = U(2:2:noDofs);

% compute stresses and displacements at mesh vertices
disp([num2str(toc),'  COMPUTING STRESSES '])

plotStressXIGA;

% plot stress

stressComp=2;
figure
clf
hold on
plot_field(node+30*[dispX dispY],elementV,'Q4',sigmaYY);
hold on
colorbar
title('Stress in x direction')
axis off
plot_mesh(node+30*[dispX dispY],elementV,'Q4','g.-');

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

disp([num2str(toc),'    SIFs COMPUTATION'])

computeInteractionIntegrals

% Compute SIFs from I integral

Knum = I.*E0/(2*(1-nu0^2)) % plain strain

% Compute the exact SIFs

Kexact = sigmato*sqrt(pi*cracklength)


%%%%%%%%%%%%

for i = 1 : size(node,1)
    x = node(i,:) ;
    [ux1,uy1] = exactDispModeI(x,E0,nu0,stressState,sigmato,xTip,seg,cracklength) ;
    ux_exact(i) = ux1;
    uy_exact(i) = uy1;
end
% ----------------------------------

% --------------------------------------------
% Plot both exact and numerical deformed shape
fac=30;
figure
hold on
h = plot(node(:,1)+fac*dispX,node(:,2)+fac*dispY,'rs');
set(h,'MarkerSize',7);
h = plot(node(:,1)+fac*ux_exact',node(:,2)+fac*uy_exact','b*');
set(h,'MarkerSize',7);
title('Exact and numerical deformed shape')
legend('XIGA','Exact')
axis equal
% --------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% post-processing for cracks

crackedMeshNURBS

