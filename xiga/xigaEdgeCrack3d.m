%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Extended Isogeometric analysis for three dimensional linear elastic
% fracture mechanics problems.
%
% Edge crack (3D) plate in tension
% Fixed at the bottom plane, uniform traction at the top plane.
%
% Vinh Phu Nguyen,
% Delft University of Technology
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

clc
clear all

global p q r controlPts weights element elementV xCr xTip 
global levelSets levelSetsB8 node

noGPs       = 2; % # of Gauss points along one direction

E0          = 1e3;  % Young modulus
nu0         = 0.3;  % Poisson ratio
sigma0      = 1;

% COMPUTE ELASTICITY MATRIX
C=zeros(6,6);
C(1:3,1:3)=E0/(1+nu0)/(1-2*nu0)*[ 1-nu0 nu0 nu0;
                                  nu0 1-nu0 nu0;
                                  nu0 nu0 1-nu0];
C(4:6,4:6)=E0/2/(1+nu0)*eye(3);

tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%  CAD input: control points and knot vectors
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

edge3dCrackC0Data;
%edge3dCrackC2Data;

% find boundary nodes for boundary conditions

fixedNodes  =  find(controlPts(:,3)==0); %bottom edge
pullNodes   =  find(controlPts(:,3)==lZ); %top edge

% Each split node is enriched by ONE function, H(x)
% Each tip node is enriched by FOUR functions, B_i(x)
% then, the total dofs is :
% total dof = numnode*nsd + numsplitnode*1*nsd + numstipnode*4*nsd
% here, two dimension, nsd = 3

noDofs = noDofs + size(split_nodes,1)*1*3 + ...
                  size(tip_nodes,1)  *4*3;

K = sparse(noDofs,noDofs);  % global stiffness matrix
u = zeros(noDofs,1);        % displacement vector
f = zeros(noDofs,1);        % external force vector

% We use fictitious nodes/control points to handle these
% additional dofs. At a H(x) enriched node, we add one fantom node and
% at tip enriched node, four fantom nodes are added. These fictitious nodes
% are numbered from the total number of true nodes, ie, from noCtrPts+1 ...

pos    = zeros(noCtrPts,1);
nsnode = 0 ;
ntnode = 0 ;

for i = 1 : noCtrPts
    if     (enrich_node(i) == 1)
        pos(i) = (noCtrPts + nsnode*1 + ntnode*4) + 1 ;
        nsnode = nsnode + 1 ;
    elseif (enrich_node(i) == 2)
        pos(i) = (noCtrPts + nsnode*1 + ntnode*4) + 1 ;
        ntnode = ntnode + 1 ;
    end
end

% essential boundary conditions

uFixed     = zeros(size(fixedNodes))';
vFixed     = zeros(size(fixedNodes))';
wFixed     = zeros(size(fixedNodes))';

udofs      = 3*fixedNodes-2;
vdofs      = 3*fixedNodes-1;
wdofs      = 3*fixedNodes;

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
    idw    = index(e,3);
    
    xiE    = elRangeU(idu,:); % [xi_i,xi_i+1]
    etaE   = elRangeV(idv,:); % [eta_j,eta_j+1]
    zetaE  = elRangeW(idw,:); % [zeta_k,zeta_k+1]
      
    % -----------------------------------------------
    % Choose Gauss quadrature rules for elements
    
    if     (ismember(e,split_elem))     % split element
        [W,Q] = quadrature(10,'GAUSS',3);
    elseif (ismember(e,tip_elem))   % tip element
        [W,Q] = quadrature(12,'GAUSS',3);
    elseif (any(intersect(tip_nodes,sctr)) ~= 0)% having tip enriched nodes
        [W,Q] = quadrature(10,'GAUSS',3);
    else
        [W,Q] = quadrature(noGPs,'GAUSS',3);
    end
    
    % Determine the position in the global matrix K
    
    sctrB = assembly3D(e,element,enrich_node,pos);
    
    % loop over Gauss points
    
    for gp=1:size(W,1)
        pt      = Q(gp,:);
      
        % compute coords in NURBS parameter space
        Xi      = parent2ParametricSpace(xiE,pt(1));
        Eta     = parent2ParametricSpace(etaE,pt(2));
        Zeta    = parent2ParametricSpace(zetaE,pt(3));
        J2      = jacobianPaPaMapping3d(xiE,etaE,zetaE);
    
        % compute derivative of basis functions w.r.t parameter coord

        [N dRdxi dRdeta dRdzeta] = NURBS3DBasisDers([Xi;Eta;Zeta],...
                                   p,q,r,uKnot,vKnot,wKnot,weights');

        % B matrix
        
        [B,J1] = BMatrixXIGA3D(e,enrich_node,N,dRdxi,dRdeta,dRdzeta,pt);
        
        % Stiffness matrix
        
        K(sctrB,sctrB) = K(sctrB,sctrB) + B'*C*B*W(gp)*J1*J2;   
    end
end

f(3*pullNodes) = 10;

disp([num2str(toc),'  APPLYING BOUNDARY  CONDITIONS'])

bcwt=mean(diag(K)); % a measure of the average  size of an element in K
% used to keep the  conditioning of the K matrix

f=f-K(:,udofs)*uFixed';  % modify the  force vector
f=f-K(:,vdofs)*vFixed';
f=f-K(:,wdofs)*wFixed';

f(udofs) = bcwt*uFixed;
f(vdofs) = bcwt*vFixed;
f(wdofs) = bcwt*wFixed;

K(udofs,:)=0;  % zero out the rows and  columns of the K matrix
K(vdofs,:)=0;
K(wdofs,:)=0;

K(:,udofs)=0;
K(:,vdofs)=0;
K(:,wdofs)=0;

K(udofs,udofs)=bcwt*speye(length(udofs));  % put ones*bcwt on the diagonal
K(vdofs,vdofs)=bcwt*speye(length(vdofs));
K(wdofs,wdofs)=bcwt*speye(length(wdofs));

% SOLVE SYSTEM

disp([num2str(toc),'  SOLVING THE SYSTEM'])

U = K\f;

% [LL UU]=lu(K);
% utemp=LL\f;
% U=UU\utemp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% POST-PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ux = U(1:3:3*noCtrPts);
Uy = U(2:3:3*noCtrPts);
Uz = U(3:3:3*noCtrPts);

% compute stresses and displacements at mesh vertices
disp([num2str(toc),'  COMPUTING STRESSES '])

vtuFile = 'edgeCrack3D';

plotStressXIGA3d

%
figure
fac = 1;
plot_mesh(node+fac*[dispX dispY dispZ],elementV,'B8','g.-',1.2);
view(3)

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,fileName,opts)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sigmaXX = zeros(size(node,1),1);
% sigmaYY = zeros(size(node,1),1);
% sigmaXY = zeros(size(node,1),1);
% 
% VTKPostProcess3d(node,elementV,'B8',...
%              [sigmaXX sigmaYY sigmaXY sigmaXX sigmaXX sigmaXX],[Ux Uy Uz]);




