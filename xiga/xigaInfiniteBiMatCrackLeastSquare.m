%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Extended Isogeometric analysis for two dimensional linear elastic
% interfacial fracture mechanics problems.
%
% Interfacial crack in an infinite plate in tension (half model)
% Dirichlet BCs = exact displacements imposed with the
% Least Square method (Luycker et al. IJNME 2011).
%
% Vinh Phu Nguyen,
% Ton Duc Thang University, Vietnam
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

global p q controlPts weights element xCrack xTips crack_node jDomainFac
global uKnot vKnot noDofs levelSets split_nodes split_elem ep CHI weakEnrFunc

nu1          = 0.3;  % Poisson ratio
E1           = 10;  % Young modulus of mat1
E2           = 1;  % of mat2
nu2          = 0.3;
stressState = 'PLANE_STRAIN';

weakEnrFunc = 1; %1=Moes func 2=absolute signed dis

vtuFile     = '../results/edgeBiMatCrackXIGA';
jDomainFac = 3;

% Elasticity matrices

Cm = elasticityMatrix(E1,nu1,stressState);
Ci = elasticityMatrix(E2,nu2,stressState);

% Constant in the 12 tip enrichment functions

G1 = E1/2/(1+nu1);       % Shear modulus for mat1
G2 = E2/2/(1+nu2);       % Shear modulus for mat2

% Kosolov constants

if strcmp(stressState,'PLANE_STRESS')
    k1 = (3-nu1)/(1+nu1);
    k2 = (3-nu2)/(1+nu2);
else
    k1 = 3-4*nu1;
    k2 = 3-4*nu2;
end

b   = (G1*(k2-1)-G2*(k1-1))/(G1*(k2+1)+G2*(k1+1));    % Second Dundur's parameter
ep  = 1/(2*pi)*log((1-b)/(1+b));                      % Material constant

cos(-ep*log(2)) - 2*ep*sin(-ep*log(2))
sin(-ep*log(2)) + 2*ep*cos(-ep*log(2))

tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%  CAD input: control points and knot vectors
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

infiniteInterfaceCrackData
%infiniteInterfaceCrackC1Data

noGPs         = p+1; % # of Gauss points along one direction
noGP1         = p+1;

% find boundary nodes for boundary conditions

bottomNodes =  find(node(:,2)==0);
rightNodes  =  find(node(:,1)==W);
topNodes    =  find(node(:,2)==W);
leftNodes   =  find(node(:,1)==0);

bottomNodesIGA =  find(controlPts(:,2)==0);
rightNodesIGA  =  find(controlPts(:,1)==W);
topNodesIGA    =  find(controlPts(:,2)==W);
leftNodesIGA   =  find(controlPts(:,1)==0);

topNodesIGA  = sort(topNodesIGA,'descend');
leftNodesIGA = sort(leftNodesIGA,'descend');

% Dirichlet control points/nodes
% should not use function "unique" because it will
% automatically sort the vector!!!

dispNodes   = [bottomNodesIGA ; rightNodesIGA(2:end);...
               topNodesIGA(2:end); leftNodesIGA(2:end-1)];
noDispNodes = length(dispNodes);

% build boundary mesh

bottomPoints    = controlPts(bottomNodesIGA,:);
rightPoints     = controlPts(rightNodesIGA,:);
topPoints       = controlPts(topNodesIGA,:);
leftPoints      = controlPts(leftNodesIGA,:);

% linear two node line elements

bottomEdgeMesh  = zeros(length(bottomNodesIGA)-1,2);
rightEdgeMesh   = zeros(length(rightNodesIGA) -1,2);
topEdgeMesh     = zeros(length(topNodesIGA)   -1,2);
leftEdgeMesh    = zeros(length(leftNodesIGA)  -1,2);

bottomEdgeMeshIGA  = zeros(noElemsU,p+1);
rightEdgeMeshIGA   = zeros(noElemsV,q+1);
topEdgeMeshIGA     = zeros(noElemsU,p+1);
leftEdgeMeshIGA    = zeros(noElemsV,q+1);

for i=1:noElemsU
    bottomEdgeMeshIGA(i,:) = bottomNodesIGA(i:i+p);
    topEdgeMeshIGA(i,:)    = topNodesIGA   (i:i+p);
end

for i=1:noElemsV
    rightEdgeMeshIGA(i,:) = rightNodesIGA(i:i+p);
    leftEdgeMeshIGA(i,:)  = leftNodesIGA (i:i+p);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Least square method for determining the control
% points displacements

disp([num2str(toc),'  LEAST SQUARE for DIRICHLET BCs'])

noBndElems = noElemsU*2 + 2*noElemsV;
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

start = max(max(bndElement));

for ie=1:noElemsU
    bndElement(ie+2*noElemsU+noElemsV,:) = [start:start+p];
    start = start + 1;
end

% final node == first node
bndElement(end,p+1) = 1;

A  = zeros(noDispNodes,noDispNodes);
bx = zeros(noDispNodes,1);
by = zeros(noDispNodes,1);

% bottom edge

noxC   = 4;

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
        
        xp    = QT*(x-xTip)';           % local coordinates
        r     = sqrt(xp(1)*xp(1)+xp(2)*xp(2));
        theta = atan2(xp(2),xp(1));
        
        % exact displacements
        
        [ux,uy] = exactDispBiMatCrack(r,theta,ep,G1,k1);
        
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
        
        xp    = QT*(x-xTip)';           % local coordinates
        r     = sqrt(xp(1)*xp(1)+xp(2)*xp(2));
        theta = atan2(xp(2),xp(1));
        
        % exact displacements
        
        [ux,uy] = exactDispBiMatCrack(r,theta,ep,G1,k1);
        
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
        
        xp    = QT*(x-xTip)';           % local coordinates
        r     = sqrt(xp(1)*xp(1)+xp(2)*xp(2));
        theta = atan2(xp(2),xp(1));
        
        % exact displacements
        
        [ux,uy] = exactDispBiMatCrack(r,theta,ep,G1,k1);
        
        bx(sctrA) = bx(sctrA) + ux*N';
        by(sctrA) = by(sctrA) + uy*N';
    end
end

% left edge

for ie=1:noElemsU
    sctr   = leftEdgeMeshIGA(ie,:);
    pts    = controlPts(sctr,:);
    sctrA  = bndElement(ie+2*noElemsU+noElemsV,:);
    xiE    = elRangeU(ie,:);
    xiArr  = linspace(xiE(1),xiE(2),noxC);
    
    for ic=1:noxC
        xi        = xiArr(ic);
        [N dNdxi] = NURBS1DBasisDers(xi,p,uKnot,weights);
        A(sctrA,sctrA) = A(sctrA,sctrA) + N'*N;
        x         = N*pts;
        
        xp    = QT*(x-xTip)';           % local coordinates
        r     = sqrt(xp(1)*xp(1)+xp(2)*xp(2));
        theta = atan2(xp(2),xp(1));
        
        % exact displacements
        
        [ux,uy] = exactDispBiMatCrack(r,theta,ep,G1,k1);
        
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

% Each split node (enrich_node(i)=1) is enriched by ONE function, H(x)
% Each tip node (enrich_node(i)=2) is enriched by FOUR functions, B_i(x)
% Each bi-mat tip node (enrich_node(i)=4) is enriched by (12+1) functions 
% Each weak discontinuity node (enrich_node(i)=3) is enriched by ONE function
% The total dofs is then:
% total dof = numnode*nsd + numsplitnode*1*nsd + numstipnode*4*nsd + ...
% here, two dimension, nsd = 2

nsd    = 2; % # of spatial dimensions
noDofs = noDofs + size(split_nodes,1)*1*nsd + ...
                  size(tip_nodes,1)*4*nsd + ... % homogeneous tip
                  size(inc_nodes,1)*1*nsd +...
                  size(itip_nodes,1)*13*nsd +... % nodes with 12 near tip & 1 weak dis.             
                  size(iTip_nodes,1)*12*nsd; % nodes with 12 near tip functions

K = sparse(noDofs,noDofs);  % global stiffness matrix
u = zeros(noDofs,1);        % displacement vector
f = zeros(noDofs,1);        % external force vector

% We use fictitious nodes/control points to handle the
% additional dofs. At a H(x) enriched node, we add one fantom node and
% at tip enriched node, four fantom nodes are added etc. 
% These fictitious nodes are numbered from the total number of true nodes, 
% ie, from numnode+1 ...

pos     = zeros(noCtrPts,1);
nsnode  = 0 ; % # of split nodes (H enriched)
ntnode  = 0 ; % # of tip enriched nodes (homogeneous crack)
ninode  = 0 ; % # of inclusion nodes
nitnode = 0 ; % # of bi-material tip/mat interface enriched
nItnode = 0 ; % # of bi-material tip enriched

for i = 1 : noCtrPts
    enrnoI = enrich_node(i);
    if     (enrnoI == 1)
        pos(i) = (noCtrPts + nsnode*1 + ntnode*4 + ninode*1 + ...
                  nitnode*13 + nItnode*12) + 1 ;
        nsnode = nsnode + 1 ;
    elseif (enrnoI == 2)
        pos(i) = (noCtrPts + nsnode*1 + ntnode*4 + ninode*1 + ...
                  nitnode*13 + nItnode*12) + 1 ;
        ntnode = ntnode + 1 ;
    elseif (enrnoI == 3)
        pos(i) = (noCtrPts + nsnode*1 + ntnode*4 + ninode*1 + ...
                  nitnode*13 + nItnode*12) + 1 ;
        ninode = ninode + 1 ;
    elseif (enrnoI == 4)
        pos(i) = (noCtrPts + nsnode*1 + ntnode*4 + ninode*1 + ...
                  nitnode*13 + nItnode*12) + 1 ;
        nitnode = nitnode + 1 ;
    elseif (enrnoI == 5)
        pos(i) = (noCtrPts + nsnode*1 + ntnode*4 + ninode*1 + ...
                  nitnode*13 + nItnode*12) + 1 ;
        nItnode = nItnode + 1 ;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Assembling system of equation
% Stiffness matrix and external force vector

noGPs         = p+1; % # of Gauss points along one direction
noGP1         = p+1;

gps = [];

disp([num2str(toc),'  ASSEMBLING THE SYSTEM'])

% Loop over elements (knot spans)

for e=1:noElems
    idu    = index(e,1);
    idv    = index(e,2);
    xiE    = elRangeU(idu,:); % [xi_i,xi_i+1]
    etaE   = elRangeV(idv,:); % [eta_j,eta_j+1]
    
    sctr   = element(e,:);          %  element scatter vector
    sctrV  = elementV(e,:);         %  element scatter vector Q4
    levelS = CHI(1,sctr);
    nn     = length(sctr);
    pts    = controlPts(sctr,:);
    
    % -----------------------------------------------
    % Choose Gauss quadrature rules for elements
    
    if     (ismember(e,split_elem))     % split element by cracks
        [W,Q] = discontQ4quad(3,levelSets(1,sctrV,1));
    elseif (ismember(e,tip_elem))   % tip element
        [W,Q] = disTipQ4quad(7,levelSets(1,sctrV),node(sctrV,:),xTip);
    elseif (any(intersect(itip_nodes,sctr)) ~= 0)% having tip enriched nodes
        [W,Q] = quadrature(10,'GAUSS',2);
    elseif (ismember(e,splitElems))     % split element by mat. interface
        [W,Q] = discontQ4quad(7,chi(1,sctrV));
    else
        [W,Q] = quadrature(noGPs,'GAUSS',2);
    end
    
    % Determine the position in the global matrix K
    
    sctrB = assembly(e,enrich_node,pos);
    
    % loop over Gauss points
    
    for gp=1:size(W,1)
        pt      = Q(gp,:);
        
        % compute coords in parameter space
        Xi      = parent2ParametricSpace(xiE, pt(1));
        Eta     = parent2ParametricSpace(etaE,pt(2));
        J2      = jacobianPaPaMapping   (xiE,etaE);
        
        % compute derivative of basis functions w.r.t parameter coord
        
        [N dRdxi dRdeta] = NURBS2DBasisDers([Xi; Eta],p,q,uKnot,vKnot,weights');
        
        % B matrix
        
        [B,J1] = BMatrixXIGA(Xi,Eta,e,enrich_node,N,dRdxi,dRdeta);
        
        % Stiffness matrix
        
        levelset = dot(N,levelS);
        x        = N * pts; % global coord of GP
        
        if levelset >= 0 
            C = Cm;
        else            
            C = Ci;
        end
        
        K(sctrB,sctrB) = K(sctrB,sctrB) + B'*C*B*W(gp)*J1*J2;
                
        gps       = [gps;x];
    end
end

disp([num2str(toc),'  APPLYING BOUNDARY  CONDITIONS'])

% SOLVE SYSTEM

disp([num2str(toc),'  SOLVING THE SYSTEM'])

applyBC

U = K\f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% POST-PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% only standard dofs
Ux = U(1:2:2*noCtrPts);
Uy = U(2:2:2*noCtrPts);

% compute stresses and displacements at mesh vertices
disp([num2str(toc),'  COMPUTING STRESSES '])

plotStressXIGAMultiMats

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,fileName,opts)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp([num2str(toc),'    SIFs COMPUTATION'])


%%%%%%%%%%%%


for i = 1 : size(node,1)
    x     = node(i,:) ;
    xp    = QT*(x-xTip)';           % local coordinates
    r     = sqrt(xp(1)*xp(1)+xp(2)*xp(2));
    theta = atan2(xp(2),xp(1));
    
    % exact displacements
    
    [ux,uy] = exactDispBiMatCrack(r,theta,ep,G1,k1);
    ux_exact(i) = ux;
    uy_exact(i) = uy;
end
% ----------------------------------

% --------------------------------------------
% Plot both exact and numerical deformed shape
dispNorm = 10/max(sqrt(Ux.^2+Uy.^2));
fac      = 1.1;
figure
hold on
plot(node(:,1)+fac*dispX,node(:,2)+fac*dispY,'r*');
%plot(node(:,1)+fac*Ux,node(:,2)+fac*Uy,'r*');
plot(node(:,1)+fac*ux_exact',node(:,2)+fac* uy_exact','b*');
%set(h,'MarkerSize',5);
title('Exact and numerical deformed shape')
legend('XIGA','Exact')
axis equal
% --------------------------------------------

figure
plot_field(node+fac*[dispX dispY],elementV,'Q4',dispY);
colorbar
figure
plot_field(node+fac*[ux_exact' uy_exact'],elementV,'Q4',uy_exact);
colorbar

stressComp=2;
figure
clf
hold on
plot_field(node+fac*[dispX dispY],elementV,'Q4',sigmaYY);
hold on
colorbar
title('Stress in x direction')
axis off
plot_mesh(node+fac*[dispX dispY],elementV,'Q4','g.-');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% post-processing for cracks

%crackedMeshNURBS

computeInteractionIntegralsBiMatCracks
