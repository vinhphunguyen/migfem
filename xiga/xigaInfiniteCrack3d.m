%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Extended Isogeometric analysis for three dimensional linear elastic
% fracture mechanics problems.
%
% Infinite crack (3D) plate.
% Exact solutions imposed on the surface boundaries with the penalty
% method.
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
addpath ../nurbs-util/
%addpath ../nurbs-geopdes/inst/

clc
clear all

global p q r controlPts weights element elementV xCr xTip 
global levelSets levelSetsB8 node


noGPs       = 3; % # of Gauss points along one direction
noGPs1      = 4;

E0          = 1e7;  % Young modulus
nu0         = 0.3;  % Poisson ratio
sigmato     = 1e4;
cracklength = 100;        % real crack length

% COMPUTE ELASTICITY MATRIX

C=zeros(6,6);
C(1:3,1:3)=E0/(1+nu0)/(1-2*nu0)*[ 1-nu0 nu0 nu0;
    nu0 1-nu0 nu0;
    nu0 nu0 1-nu0];
C(4:6,4:6)=E0/(1+nu0)*eye(3);

tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%  CAD input: control points and knot vectors
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp([num2str(toc),'  INPUT and ENRICHMENT '])

%infinite3dCrackC0Data;
infinite3dCrackC1Data;
%infinite3dCrackC2Data;

% find boundary nodes for boundary conditions

bottomNodes =  find(controlPts(:,3)==0);
rightNodes  =  find(controlPts(:,1)==lX);
topNodes    =  find(controlPts(:,3)==lZ);
leftNodes   =  find(controlPts(:,1)==0);

% generate surface mesh from the above node sets

[botElems,botIndex]     = surfaceMesh (uKnot,vKnot,bottomNodes,p,q,...
    noPtsX,noPtsY,elRangeU,elRangeV,elConnU,elConnV);

[topElems,topIndex]     = surfaceMesh (uKnot,vKnot,topNodes,p,q,...
    noPtsX,noPtsY,elRangeU,elRangeV,elConnU,elConnV);

[rightElems,rightIndex] = surfaceMesh (vKnot,wKnot,rightNodes,q,r,...
    noPtsY,noPtsZ,elRangeV,elRangeW,elConnV,elConnW);

[leftElems,leftIndex]   = surfaceMesh (vKnot,wKnot,leftNodes,q,r,...
    noPtsY,noPtsZ,elRangeV,elRangeW,elConnV,elConnW);

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
    
    if     (ismember(e,split_elem))             % split element
        [W,Q] = quadrature(10,'GAUSS',3);
    elseif (ismember(e,tip_elem))               % tip element
        [W,Q] = quadrature(12,'GAUSS',3);
    elseif (any(intersect(tip_nodes,sctr)) ~= 0)% having tip enriched nodes
        [W,Q] = quadrature(12,'GAUSS',3);
    else
        [W,Q] = quadrature(noGPs,'GAUSS',3);
    end
    
    % Determine the position in the global matrix K
    
    sctrB = assembly3D(e,element,enrich_node,pos);
    
    % loop over Gauss points
    
    for gp=1:size(W,1)
        pt      = Q(gp,:);
        
        % compute coords in NURBS parameter space
        
        Xi      = parent2ParametricSpace(xiE,  pt(1));
        Eta     = parent2ParametricSpace(etaE, pt(2));
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

disp([num2str(toc),'  APPLYING BOUNDARY  CONDITIONS'])

% If penalty method is use, then one must modify the stiffness matrix and
% the nodal force vector
% $K_{ij}$ = Kij - alpha \int phi_i phi_j d \gamma_u
% fj       = fj  - alpha \int phi_i u_bar d \gamma_u

fu = zeros(noDofs,1);
k  = zeros(noDofs,noDofs);

[W1,Q1] = quadrature(noGPs1, 'GAUSS', 2);

S   = [1 0 0; 0 1 0; 0 0 1];

% Loop over elements of the bottom surface

for e=1:length(botElems)
    idu    = botIndex(e,1);
    idv    = botIndex(e,2);
    xiE    = elRangeU(idu,:); % [xi_i,xi_i+1]
    etaE   = elRangeV(idv,:); % [eta_j,eta_j+1]
    sctr   = botElems(e,:);
    pts    = controlPts(sctr,:);
    pts    = pts(:,1:2);            
    
    le     = length(sctr);
    en     = zeros(1,3*le);
    force  = zeros(1,3*le);
    
    Phi        = levelSets(sctr,:);
    normalPhi  = Phi(:,1);
    tangentPhi = Phi(:,2);
        
    % loop over Gauss points
    
    for gp=1:size(W1,1)
        pt      = Q1(gp,:);
        wt      = W1(gp);
        Xi      = parent2ParametricSpace(xiE, pt(1));
        Eta     = parent2ParametricSpace(etaE,pt(2));
        J2      = jacobianPaPaMapping(xiE,etaE);
        
        [N dRdxi dRdeta] = NURBS2DBasisDers([Xi;Eta],p,q,uKnot,vKnot,weights');
        
        % compute the jacobian of physical and parameter domain mapping
                            
        jacob   = pts'*[dRdxi' dRdeta'];
        J1      = det(jacob);
        
        % compute exact displacements
       
        phi   = N * normalPhi;
        psi   = N * tangentPhi;

        [ux,uy,uz] = exactGriffith3D(phi,psi,E0,nu0,sigmato,cracklength);
                
        for j = 1 : le
            tem = 3*sctr(j);
            j3  = 3*j;
            Nj  = N(j);
            en(j3-2) = tem-2;
            en(j3-1) = tem-1;
            en(j3  ) = tem;
            
            force(j3-2) = S(1,1)*Nj*ux;
            force(j3-1) = S(2,2)*Nj*uy;
            force(j3  ) = S(3,3)*Nj*uz;
        end
        
        fac    = J1 * J2 * wt;
        
        fu(en) = fu(en) + fac * force';
        
        for i = 1:le
            id=[3*sctr(i)-2 3*sctr(i)-1 3*sctr(i)];
            for j=1:le
                jd=[3*sctr(j)-2 3*sctr(j)-1 3*sctr(j)];
                k(id,jd) = k(id,jd) + fac * N(i)*N(j)*S ;
            end
        end
    end % end of loop over Gauss points
end % end of loop over boundary elements

% Loop over elements of the top surface

for e=1:length(topElems)
    idu    = topIndex(e,1);
    idv    = topIndex(e,2);
    xiE    = elRangeU(idu,:); % [xi_i,xi_i+1]
    etaE   = elRangeV(idv,:); % [eta_j,eta_j+1]
    sctr   = topElems(e,:);
    pts    = controlPts(sctr,:);
    pts    = pts(:,1:2);          % get rid of z-coord  
                
    le    = length(sctr);
    en    = zeros(1,3*le);
    force = zeros(1,3*le);
    
    Phi        = levelSets(sctr,:);
    normalPhi  = Phi(:,1);
    tangentPhi = Phi(:,2);
    
    % loop over Gauss points
    
    for gp=1:size(W1,1)
        pt      = Q1(gp,:);
        wt      = W1(gp);
        Xi      = parent2ParametricSpace(xiE, pt(1));
        Eta     = parent2ParametricSpace(etaE,pt(2));
        J2      = jacobianPaPaMapping(xiE,etaE);
        
        [N dRdxi dRdeta] = NURBS2DBasisDers([Xi;Eta],p,q,uKnot,vKnot,weights');
        
        % compute the jacobian of physical and parameter domain mapping
        
        x       = N   * pts; % global coord of GP              
        jacob   = pts'*[dRdxi' dRdeta'];
        J1      = det(jacob);
        
        % compute exact displacements
        
        phi   = N * normalPhi;
        psi   = N * tangentPhi;

        [ux,uy,uz] = exactGriffith3D(phi,psi,E0,nu0,sigmato,cracklength);
                
        for j = 1 : le
            tem = 3*sctr(j);
            j3  = 3*j;
            Nj  = N(j);
            en(j3-2) = tem-2;
            en(j3-1) = tem-1;
            en(j3  ) = tem;
            
            force(j3-2) = S(1,1)*Nj*ux;
            force(j3-1) = S(2,2)*Nj*uy;
            force(j3  ) = S(3,3)*Nj*uz;
        end
        
        fac    = J1 * J2 * wt;
        fu(en) = fu(en) + fac * force';
        
        for i = 1:le
            id=[3*sctr(i)-2 3*sctr(i)-1 3*sctr(i)];
            for j=1:le
                jd=[3*sctr(j)-2 3*sctr(j)-1 3*sctr(j)];
                k(id,jd) = k(id,jd) + fac * N(i)*N(j)*S ;
            end
        end
    end % end of loop over Gauss points
end % end of loop over boundary elements

% Loop over elements of the right surface

for e=1:length(rightElems)
    idu    = rightIndex(e,1);
    idv    = rightIndex(e,2);
    xiE    = elRangeV(idu,:); % [xi_i,xi_i+1]
    etaE   = elRangeW(idv,:); % [eta_j,eta_j+1]
    sctr   = rightElems(e,:);
    pts    = controlPts(sctr,:);
    pts    = pts(:,2:3);     % get rid of x-coord       
       
    le    = length(sctr);
    en    = zeros(1,3*le);
    force = zeros(1,3*le);
        
    Phi        = levelSets(sctr,:);
    normalPhi  = Phi(:,1);
    tangentPhi = Phi(:,2);
    
    % loop over Gauss points
    
    for gp=1:size(W1,1)
        pt      = Q1(gp,:);
        wt      = W1(gp);
        Xi      = parent2ParametricSpace(xiE, pt(1));
        Eta     = parent2ParametricSpace(etaE,pt(2));
        J2      = jacobianPaPaMapping(xiE,etaE);
        
        [N dRdxi dRdeta] = NURBS2DBasisDers([Xi;Eta],q,r,vKnot,wKnot,weights');
        
        % compute the jacobian of physical and parameter domain mapping
        
        x       = N   * pts; % global coord of GP              
        jacob   = pts'*[dRdxi' dRdeta'];
        J1      = det(jacob);
        
        % compute exact displacements
        
        phi   = N * normalPhi;
        psi   = N * tangentPhi;

        [ux,uy,uz] = exactGriffith3D(phi,psi,E0,nu0,sigmato,cracklength);
                
        for j = 1 : le
            tem = 3*sctr(j);
            j3  = 3*j;
            Nj  = N(j);
            en(j3-2) = tem-2;
            en(j3-1) = tem-1;
            en(j3  ) = tem;
            
            force(j3-2) = S(1,1)*Nj*ux;
            force(j3-1) = S(2,2)*Nj*uy;
            force(j3  ) = S(3,3)*Nj*uz;
        end
        
        fac    = J1 * J2 * wt;        
        fu(en) = fu(en) + fac * force';
        
        for i = 1:le
            id=[3*sctr(i)-2 3*sctr(i)-1 3*sctr(i)];
            for j=1:le
                jd=[3*sctr(j)-2 3*sctr(j)-1 3*sctr(j)];
                k(id,jd) = k(id,jd) + fac * N(i)*N(j)*S ;
            end
        end
    end % end of loop over Gauss points
end % end of loop over boundary elements

% External force on the left surface

% normal vector to this surface
nx = -1;
ny = 0;
nz = 0;

% for e=1:length(leftElems)
%     idu    = leftIndex(e,1);
%     idv    = leftIndex(e,2);
%     xiE    = elRangeV(idu,:); % [xi_i,xi_i+1]
%     etaE   = elRangeW(idv,:); % [eta_j,eta_j+1]
%     sctr   = leftElems(e,:);
%     pts    = controlPts(sctr,:);
%     pts    = pts(:,2:3);            
%     
%     le     = length(sctr);
%     en     = zeros(1,3*le);
%     force  = zeros(1,3*le);
%     
%     Phi        = levelSets(sctr,:);
%     normalPhi  = Phi(:,1);
%     tangentPhi = Phi(:,2);
%     
%     sctrx  = 3*sctr-2;
%     sctry  = 3*sctr-1;
%     sctrz  = 3*sctr;
% 
%     % loop over Gauss points
%     
%     for gp=1:size(W1,1)
%         pt      = Q1(gp,:);
%         wt      = W1(gp);
%         Xi      = parent2ParametricSpace(xiE, pt(1));
%         Eta     = parent2ParametricSpace(etaE,pt(2));
%         J2      = jacobianPaPaMapping(xiE,etaE);
%         
%         [N dRdxi dRdeta] = NURBS2DBasisDers([Xi;Eta],q,r,vKnot,wKnot,weights');
%         
%         % compute the jacobian of physical and parameter domain mapping
%                             
%         jacob   = pts'*[dRdxi' dRdeta'];
%         J1      = det(jacob);
%         
%         % compute exact displacements
%        
%         phi   = N * normalPhi;
%         psi   = N * tangentPhi;
%                 
%         [sigmaxx,sigmayy,sigmazz,sigmayz,sigmaxz,sigmaxy] = ...
%             exactGriffithStress3D(phi,psi,E0,nu0,sigmato,cracklength);
%         
%         tx = sigmaxx*nx + sigmaxy*ny + sigmaxz*nz;
%         ty = sigmaxy*nx + sigmayy*ny + sigmayz*nz;
%         tz = sigmaxz*nx + sigmayz*ny + sigmazz*nz;
%        
%         fac = J1 * J2 * wt;
%         
%         f(sctrx) = f(sctrx) + N' * tx * fac;
%         f(sctry) = f(sctry) + N' * ty * fac;
%         f(sctrz) = f(sctrz) + N' * tz * fac;
%     end % end of loop over Gauss points
% end % end of loop over boundary elements


% modified system of equations

alpha = 0.05*max(max(C)) ;           % penalty number
alpha = 1e10;
f     = f - alpha*fu;
K     = K - alpha*k;

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

vtuFile = 'infiniteCrack3D';

plotStressXIGA3d

%
figure
fac = 30;
plot_mesh(node+fac*[dispX dispY dispZ],elementV,'B8','g.-');
view(3)

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,fileName,opts)


for i = 1 : size(node,1) 
    phi = levelSetsB8(i,1);
    psi = levelSetsB8(i,2);
    [ux1,uy1,uz1] = exactGriffith3D(phi,psi,E0,nu0,sigmato,cracklength) ;
    ux_exact(i) = ux1;
    uy_exact(i) = uy1;
    uz_exact(i) = uz1;
end
% ----------------------------------

% --------------------------------------------
% Plot both exact and numerical deformed shape
fac=30;
figure
hold on
h = plot3(node(:,1)+fac*dispX,...
          node(:,2)+fac*dispY,...
          node(:,3)+fac*dispZ,'rs');
set(h,'MarkerSize',7);
h = plot3(node(:,1)+fac*ux_exact',...
          node(:,2)+fac*uy_exact',...
          node(:,3)+fac*uz_exact','b*');
set(h,'MarkerSize',7);
title('Exact and numerical deformed shape')
legend('XIGA','Exact')
set(gcf, 'color', 'white');
axis equal
axis off
view(3)
% --------------------------------------------

vtuFile = 'infiniteCrack3dExact';
VTKPostProcess3d(node,elementV,'B8',vtuFile,...
    [sigmaXX sigmaYY sigmaZZ sigmaXY sigmaYZ sigmaZX], ...
    [ux_exact' uy_exact' uz_exact'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


