%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Isogeometric analysis for Kirchoff plate problems.
%
% Rotation-free thin plates. Fully clamped or simply supported 
% circular plates.
%
% Vinh Phu Nguyen,
% Cardiff University
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

tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%  CAD input: control points and knot vectors
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

circularPlateData

% constitutive matrix

D  = E*t^3/(12*(1-nu^2));
C  = D*[1 nu 0;nu 1 0; 0 0 0.5*(1-nu)];

% build connectivity ...

generateIGA2DMesh

noCtrPts       = noPtsX   * noPtsY;
noDofs         = noCtrPts * 1;

% find boundary nodes for boundary conditions

EPS = 1e-8;
bottomNodes  =  1:noPtsX;
topNodes     =  noPtsX*(noPtsY-1)+1:noCtrPts;
leftNodes    =  1:noPtsX:noPtsX*(noPtsY-1)+1;
rightNodes   =  noPtsX:noPtsX:noCtrPts;

fixedNodes  =  unique([bottomNodes;rightNodes;topNodes;leftNodes]);

if clamped
nextToBotNodes = noPtsX+2:2*noPtsX-1;
nextToRgtNodes = 2*noPtsX-1:noPtsX:noPtsX*(noPtsY-1)-1;
nextToTopNodes = noPtsX*(noPtsY-2)+2:noPtsX*(noPtsY-1)-1;
nextToLefNodes = noPtsX+2:noPtsX:noPtsX*(noPtsY-2)+2;

nextNodes      = unique([nextToBotNodes';nextToRgtNodes';...
                         nextToTopNodes';nextToLefNodes']);

fixedNodes     = [fixedNodes; nextNodes(:)];
end

plot(controlPts(fixedNodes,1),controlPts(fixedNodes,2),...
    'bs','MarkerEdgeColor','r','MarkerSize',14);


% initialization

K = sparse(noDofs,noDofs);  % global stiffness matrix
f = zeros(noDofs,1);        % external force vector

% essential boundary conditions

uFixed     = zeros(size(fixedNodes))';
udofs      = fixedNodes;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gauss quadrature rule
noGPs = p+1;
noGpEle = noGPs^2;
[W,Q]=quadrature(  noGPs, 'GAUSS', 2 ); % noGPs x noGPs point quadrature

% Assembling system of equation
% Stiffness matrix and external force vector

disp([num2str(toc),'  ASSEMBLING THE SYSTEM']);

% Loop over elements (knot spans)

for e=1:noElems
    idu    = index(e,1);
    idv    = index(e,2);
    xiE    = elRangeU(idu,:); % [xi_i,xi_i+1]
    etaE   = elRangeV(idv,:); % [eta_j,eta_j+1]
    
    sctr   = element(e,:);          %  element scatter vector
    nn     = length(sctr);
    pts    = controlPts(sctr,:);
    
    % loop over Gauss points
    
    for gp=1:size(W,1)
        pt      = Q(gp,:);
        wt      = W(gp);
        
        % compute coords in parameter space
        Xi      = parent2ParametricSpace(xiE,pt(1));
        Eta     = parent2ParametricSpace(etaE,pt(2));
        J2      = jacobianPaPaMapping(xiE,etaE);
        
        % shape functions, first and second derivatives w.r.t natural coords
        
        [R dRdxi dRdeta dR2dxi dR2det dR2dxe] = ...
            NURBS2DBasis2ndDers([Xi; Eta],p,q,uKnot,vKnot,weights');
        
        % compute the jacobian of physical and parameter domain mapping
        % then the derivative w.r.t spatial physical coordinates
        
        jacob  = [dRdxi; dRdeta]          * pts; % 2x2 matrix
        jacob2 = [dR2dxi; dR2det; dR2dxe] * pts; % 3x2 matrix
        
        J1    = det(jacob);
        
        dxdxi = jacob(1,1); dydxi = jacob(1,2);
        dxdet = jacob(2,1); dydet = jacob(2,2);
        
        j33   = [dxdxi^2     dydxi^2     2*dxdxi*dydxi;
            dxdet^2     dydet^2     2*dxdet*dydet;
            dxdxi*dxdet dydxi*dydet dxdxi*dydet+dxdet*dydxi];
        
        % Jacobian inverse and spatial 1st and 2nd derivatives
        
        invJacob   = inv(jacob);
        dRdx       = invJacob*[dRdxi;dRdeta];
        dR2dx      = inv(j33)*([dR2dxi; dR2det; dR2dxe]-jacob2*dRdx);
        
        % B matrix
        
        B          = dR2dx;
        B(3,:)     = B(3,:)*2;
        
        % compute elementary stiffness matrix and
        % assemble it to the global matrix
        
        K(sctr,sctr) = K(sctr,sctr) + B' * C * B * J1 * J2 * wt;
        f(sctr)      = f(sctr)      + q0 * R' * J1 * J2 * wt;
    end
end

disp([num2str(toc),'  APPLYING BOUNDARY  CONDITIONS']);
bcwt=mean(diag(K)); % a measure of the average  size of an element in K
% used to keep the  conditioning of the K matrix

f=f-K(:,udofs)*uFixed';  % modify the  force vector
f(udofs) = bcwt*uFixed;
K(udofs,:)=0;  % zero out the rows and  columns of the K matrix
K(:,udofs)=0;
K(udofs,udofs)=bcwt*speye(length(udofs));  % put ones*bcwt on the diagonal

% SOLVE SYSTEM
disp([num2str(toc),'  SOLVING THE SYSTEM']);
U=K\f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% POST-PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Post-processing using triangulation of GPs

% xx    = zeros(noGpEle*noElems,2);  % global coords of Gauss points
% u     = zeros(noGpEle*noElems,1);  % global coords of Gauss points
%
%
% id=1;
% for e=1:noElems
%     sctr   = element(e,:);         %  element scatter vector
%     % loop over Gauss points
%     pts    = controlPts(sctr,:);
%     for gp=1:size(W,1)
%         pt      = Q(gp,:);
%         wt      = W(gp);
%
%         [R dRdxi dRdeta] = ...
%           NURBS2DBasisDers([Xi; Eta],p,q,uKnot,vKnot,weights');
%
%         xx(id,:)    = R  * pts;
%         u(id,:)     = R  * U(sctr);
%
%        id = id + 1;
%     end
% end
%
% figure
% plot3(xx(:,1),xx(:,2),u,'.');
%
% tri   = delaunay(xx);
% noTri = size(tri,1);
%
% X = zeros(3,noTri);
% Y = zeros(3,noTri);
% Z = zeros(3,noTri);
%
% for i = 1:noTri
%     sctr   = tri(i,:);
%     X(:,i) = xx(sctr,1);
%     Y(:,i) = xx(sctr,2);
%     Z(:,i) = u   (sctr,1);
% end
%
% figure
% fill3(X,Y,Z,Z);
% colorbar
% title('deflection')
% axis on

%% Visualization using a Q4 visualization mesh

buildVisualizationMesh;
%plot_mesh(node,elementV,'Q4','g.-',1);

disp   = zeros(noElems,size(elementV,2));

for e=1:noElems
    idu    = index(e,1);
    idv    = index(e,2);
    xiE    = elRangeU(idu,:); % [xi_i,xi_i+1]
    etaE   = elRangeV(idv,:); % [eta_j,eta_j+1]
    
    sctr   = element(e,:);         %  element scatter vector
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
            
            disp(e,gp)    = N * U(sctr);
            
            gp = gp +1;
        end
    end
end

X = zeros(4,noElemsV);
Y = zeros(4,noElemsV);
Z = disp';

for i = 1:size(elementV,1)
    sctr   = elementV(i,:);
    X(:,i) = node(sctr,1);
    Y(:,i) = node(sctr,2);
end

figure
fill3(X,Y,Z,Z);
colorbar
title('deflection')
axis off


opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,fileName,opts)


% Exact solution check

wbar = min(Z(:))
wext = q0*r^4/64/D  % CCCC
%wext = 4.06235; % SSSS







