%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Extended IGA for two dimensional elasticity problems.
%
% Infinite plate with a circular inclusion (Tri's Msc thesis).
% The enrichment function is the modified abs function proposed
% by Moes et al 2003 or simply the absolute signed distance function (Sukumar).
%
% Vinh Phu Nguyen, March 2012
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


E1           = 1000;  % Young modulus of matrix
nu1          = 0.3;   % Poisson ratio
E2           = 1;     % of inclusion
nu2          = 0.3;
stressState = 'PLANE_STRAIN';
sigma0      = 10;

levelMax    = 3;

weakEnrFunc = 2;

vtuFile     = '../results/infinitePlateInclusionXIGA';

% Elasticity matrix

Cm = elasticityMatrix(E1,nu1,stressState);
Ci = elasticityMatrix(E2,nu2,stressState);

tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%  CAD input: control points and knot vectors
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plateHoleData

noGPs          = p+1; % # of GPs, assume p>=q

% Boundary nodes
% 1. Displacement conditions (symmetry):
%    - On bottom edge (y=0) with u_y = 0
%    - On left  edge (x=0) with u_x = 0

% 2. Traction conditions
%    - On right edge (x=L) with t = (sigma_x,sigma_xy)
%    - On top edge  (y=L) with t = (sigma_xy,sigma_y)

% find boundary nodes for boundary conditions

bottomNodes = find(controlPts(:,2)==0)';
rightNodes  = find(controlPts(:,1)==L)';
leftNodes   = find(controlPts(:,1)==0)';
topNodes    = find(controlPts(:,2)==2*L)';

% essential boundary conditions

uNodes = bottomNodes;
vNodes = [bottomNodes];

uFixed     = zeros(size(uNodes));
vFixed     = zeros(size(bottomNodes));

udofs = 2*uNodes-1; % global indecies  of the fixed x disps
vdofs = 2*vNodes;   % global indecies  of the fixed y disps

% build a 1D meshes for the right and top edges over which
% a traction is applied

rightPoints   = controlPts(rightNodes,:);
leftPoints    = controlPts(leftNodes,:);
topPoints     = controlPts(topNodes,:);
bottomPoints  = controlPts(bottomNodes,:);
rightEdgeMesh = zeros(noElemsV,q+1);
leftEdgeMesh  = zeros(noElemsV,q+1);
topEdgeMesh   = zeros(noElemsU,p+1);

for i=1:noElemsV
    rightEdgeMesh(i,:) = rightNodes(i:i+q);
    leftEdgeMesh(i,:)  = leftNodes(i:i+q);
end

for i=1:noElemsU
    topEdgeMesh(i,:)    = topNodes(i:i+p);
    bottomEdgeMesh(i,:) = bottomNodes(i:i+p);
end

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
    
    W      = [];
    Q      = [];
    
    % -----------------------------------------------
    % Choose Gauss quadrature rules for elements
    
    voidId     = find(splitElems(:,1)==e);
    if (isempty(voidId) == 0)     % split element
        [aa] = hierarchicalGaussQuad(noGPs,CHI(elementV(e,:)),...
            coord,node(elementV(e,:),:),splitElems(voidId,2),0);
        %[W,Q] = discontQ4quad(7,ls(1,elementV(e,:),1));
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
        
        [R, dRdxi, dRdeta] = NURBS2DBasisDers([Xi;Eta],p,q,uKnot,vKnot,weights');
        
        % compute the jacobian of physical and parameter domain mapping
        % then the derivative w.r.t spatial physical coordinates
        % Jacobian matrix
        
        x         = R * pts; % global coord of GP
        gps       = [gps;x];
        
        xc       = VOID(iVoid,1);
        yc       = VOID(iVoid,2);
        rc       = VOID(iVoid,3);
        levelset = sqrt((x(1)-xc)^2+(x(2)-yc)^2)-rc;
        
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

% Computing external force
%%
[W1,Q1] = quadrature(noGPs, 'GAUSS', 1 );

% loop over top edge

for e=1:noElemsU
    xiE   = elRangeU(e,:); % [xi_i,xi_i+1]
    conn  = elConnU(e,:);
    pts   = topPoints(conn,:);
    sctrx = topEdgeMesh(e,:);
    sctry = 2*sctrx;
    
    % loop over Gauss points
    for gp=1:size(W1,1)
        xi      = Q1(gp,:);
        wt      = W1(gp);
        Xi      = 0.5 * ( ( xiE(2) - xiE(1) ) * xi + xiE(2) + xiE(1));
        J2      = 0.5 * ( xiE(2) - xiE(1) );
        
        [N, dNdxi] = NURBS1DBasisDers(Xi,p,uKnot,weights);
        
        x        = N     * pts; % global coord of GP
        jacob1   = dNdxi * pts;
        J11      = norm (jacob1);
        
        tx       = 0;
        ty       = sigma0;
        fac      = J11 * J2 * wt;
        f(sctrx) = f(sctrx) + N' * tx * fac;
        f(sctry) = f(sctry) + N' * ty * fac;
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

Ux    = U(1:noCtrPts);
Uy    = U(noCtrPts+1:noDofs);

%%

% Plot mesh and enriched nodes to check
figure
hold on
plot_mesh(node,elementV,'Q4','b-',1.2);
% plot the circle
theta = 0:0.01:2*pi;
xo = xc + r*cos(theta) ;
yo = yc + r*sin(theta) ;
plot(xo,yo,'k-','Linewidth',1.9);
% plot elements cut by the circle
plot_mesh(node,elementV(splitElems(:,1),:),'Q4','r-',1.2);
plot(gps(:,1),gps(:,2),'+');
n1 = plot(controlPts(inc_nodes,1),controlPts(inc_nodes,2),'r*');
set(n1,'MarkerSize',16,'LineWidth',1.07);

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
% exportfig(gcf,fileName,opts)

%%%%%%%%%%%%%%%%%%%%%
%plot stress field

stress = zeros(noElems,size(elementV,2),4);
disp   = zeros(noElems,size(elementV,2),2);

for e=1:noElems
    idu    = index(e,1);
    idv    = index(e,2);
    xiE    = elRangeU(idu,:); % [xi_i,xi_i+1]
    etaE   = elRangeV(idv,:); % [eta_j,eta_j+1]
    
    sctr   = element(e,:);          %  element scatter vector
    sctrV  = elementV(e,:);        %  element scatter vector
    nn     = length(sctr);
    
    pts    = controlPts(sctr,:);
    
    uspan = FindSpan(noPtsX-1,p,xiE(1),uKnot);
    vspan = FindSpan(noPtsY-1,q,etaE(1),vKnot);
    
    elemDisp  = element_disp(e,pos,enrich_node,U);
    
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
            
            [B,w] = BMatrixXIGA(Xi,Eta,e,enrich_node,N,dRdxi,dRdeta);
            [exN] = NMatrixXIGA(e,enrich_node,N);
            
            if (chi(sctrV(gp)) >= 0)
                C = Cm;
            else
                C = Ci;
            end
            
            strain          = B*elemDisp;
            sigma           = C*strain;
            stress(e,gp,1:3)= sigma;
            
            % von Mises stress
            stress(e,gp,4)  = sqrt(sigma(1)^2+sigma(2)^2-...
                                   sigma(1)*sigma(2)+3*sigma(3)^2);
            
            % the following is incorrect since
            % enriched dofs are not included!!!
            %disp(e,gp,:)    = N*[Ux(sctr) Uy(sctr)];
            
            disp(e,gp,:)    = exN*[elemDisp(1:2:end) ...
                                   elemDisp(2:2:end)];
            gp = gp +1;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% export to VTK format to plot in Mayavi or Paraview

nNodes = size(node,1);

Sxx = zeros(nNodes,2); Sxy = Sxx; Syy = Sxx; Svm = Sxx;

dispX = zeros(size(node,1),1);
dispY = zeros(size(node,1),1);

for e=1:size(elementV,1)
    connect = elementV(e,:);
    for in=1:4
        nNode = connect(in);
        
        Sxx(nNode,:) = Sxx(nNode,:) + [stress(e,in,1) 1];
        Sxy(nNode,:) = Sxy(nNode,:) + [stress(e,in,3) 1];
        Syy(nNode,:) = Syy(nNode,:) + [stress(e,in,2) 1];
        Svm(nNode,:) = Svm(nNode,:) + [stress(e,in,4) 1];
        
        dispX(nNode) = disp(e,in,1);
        dispY(nNode) = disp(e,in,2);
    end
end

% Average nodal stress values
Sxx(:,1) = Sxx(:,1)./Sxx(:,2); Sxx(:,2) = [];
Sxy(:,1) = Sxy(:,1)./Sxy(:,2); Sxy(:,2) = [];
Syy(:,1) = Syy(:,1)./Syy(:,2); Syy(:,2) = [];
Svm(:,1) = Svm(:,1)./Svm(:,2); Svm(:,2) = [];


VTKPostProcess(node,elementV,2,'Quad4',vtuFile,...
    [Sxx(:,1) Syy(:,1) Sxy(:,1) Svm(:,1)],[dispX dispY]);

fac =5;

stressComp=4;
% figure
% clf
% plot_field(node,elementV,'Q4',Sxx(:,1));
% hold on
% colorbar
% title('von Mises stress')
% axis off
%plot_mesh(node,elementV,'Q4','g.-');

figure
clf
plot_field(node+fac*[dispX dispY],elementV,'Q4',disp(:,:,1));
hold on
colorbar
title('Displacement in x direction')
axis off

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,fileName,opts)








