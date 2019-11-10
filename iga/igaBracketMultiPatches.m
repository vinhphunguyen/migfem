%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Isogeometric analysis for two dimensional elasticity problems.
% Illustration of a multi-patch IGA code.
% Bracket
%
% Vinh Phu Nguyen, March 2012
% Delft University of Technology
% nvinhphu@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('~/code/xfem-efg-matlab/fem_util');
addpath C_files/

clc
clear all

E           = 1e5;  % Young modulus
nu          = 0.3;  % Poisson’s ratio
stressState = 'PLANE_STRESS';
a           = 0.5;

% Elasticity matrix

if ( strcmp(stressState,'PLANE_STRESS') )
    C=E/(1-nu^2)*[ 1      nu          0;
        nu     1          0 ;
        0     0  0.5*(1-nu) ];
else
    C=E/(1+nu)/(1-2*nu)*[ 1-nu  nu     0;
        nu    1-nu   0;
        0     0  0.5-nu ];
end

tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%  CAD input: control points and knot vectors
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bracketData

noGPs  = 3;
noGPs1 = noGPs;

noCtrPts       = max(max(patches(noPatches).element));
noDofs         = noCtrPts * 2;

% find boundary nodes for boundary conditions

leftNodes     = nodePattern2(end,:);
bottomNodes   = nodePattern3(end,:);
topNodes      = nodePattern (end,:);
rightNodes    = nodePattern4(end,:);

idx           = nodePattern(end,:);

leftPoints    = patch2.controlPts(idx,:);
topPoints     = patch1.controlPts(idx,:);
rightPoints   = patch4.controlPts(idx,:);
bottomPoints  = patch3.controlPts(idx,:);

leftEdgeMesh   = zeros(patch1.noElemsV,patch1.q+1);
rightEdgeMesh  = zeros(patch1.noElemsV,patch1.q+1);
topEdgeMesh    = zeros(patch1.noElemsV,patch1.q+1);
bottomEdgeMesh = zeros(patch1.noElemsV,patch1.q+1);

for i=1:patch1.noElemsV
    leftEdgeMesh(i,:)   = leftNodes  (i:i+patch1.q);
    rightEdgeMesh(i,:)  = rightNodes (i:i+patch1.q);
    topEdgeMesh(i,:)    = topNodes   (i:i+patch1.q);
    bottomEdgeMesh(i,:) = bottomNodes(i:i+patch1.q);
end

% initialization

K = sparse(noDofs,noDofs);  % global stiffness matrix
u = zeros(noDofs,1);        % displacement vector
f = zeros(noDofs,1);        % external force vector

% essential boundary conditions, exact displacement is used
% as the left edge
vNodes = [leftNodes(end) rightNodes(1)];
uFixed     = [zeros(size(leftNodes)) ones(size(rightNodes))];
vFixed     = zeros(size(vNodes));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Gauss quadrature rule
[W,Q]=quadrature(  noGPs, 'GAUSS', 2 ); % noGPs x noGPs point quadrature

% Assembling system of equation
% Stiffness matrix and external force vector

disp([num2str(toc),'  ASSEMBLING THE SYSTEM'])

% Loop over patches

for ip=1:noPatches
    index      = patches(ip).index;
    elRangeU   = patches(ip).elRangeU;
    elRangeV   = patches(ip).elRangeV;
    element    = patches(ip).element;    
    elementL   = patches(ip).elementLocal;    
    uKnot      = patches(ip).uKnot;
    vKnot      = patches(ip).vKnot;
    controlPts = patches(ip).controlPts;
    weights    = patches(ip).weights;
    noElems    = patches(ip).noElemsU * patches(ip).noElemsV;
    p          = patches(ip).p;
    q          = patches(ip).q;
    
    % Loop over elements (knot spans)
    
    for e=1:noElems
        idu    = index(e,1);
        idv    = index(e,2);
        xiE    = elRangeU(idu,:); % [xi_i,xi_i+1]
        etaE   = elRangeV(idv,:); % [eta_j,eta_j+1]
        
        sctr   = element(e,:);         %  element scatter vector
        sctrL  = elementL(e,:);         %  element scatter vector
        sctrB  = [sctr sctr+noCtrPts]; %  vector that scatters a B matrix
        nn     = length(sctr);
        pts    = controlPts(sctrL,:);
        
        B      = zeros(3,2*nn);
        
        % loop over Gauss points
        
        for gp=1:size(W,1)
            pt      = Q(gp,:);
            wt      = W(gp);
            
            % compute coords in parameter space
            Xi      = parent2ParametricSpace(xiE, pt(1));
            Eta     = parent2ParametricSpace(etaE,pt(2));
            J2      = jacobianPaPaMapping(xiE,etaE);
            
            % compute derivatives of shape functions
            [R dRdxi dRdeta] = NURBS2DBasisDers([Xi;Eta],p,q,uKnot,vKnot,weights');
            
            % compute the jacobian of physical and parameter domain mapping
            % then the derivative w.r.t spatial physical coordinates
            
            jacob      = pts' * [dRdxi' dRdeta'];
            J1         = det(jacob);                                    
            invJacob   = inv(jacob);
            dRdx       = [dRdxi' dRdeta'] * invJacob;
            
            % B matrix

            B(1,1:nn)       = dRdx(:,1)';
            B(2,nn+1:2*nn)  = dRdx(:,2)';
            B(3,1:nn)       = dRdx(:,2)';
            B(3,nn+1:2*nn)  = dRdx(:,1)';
            
            % elementary stiffness matrix and assemble it to the global matrix
            
            K(sctrB,sctrB) = K(sctrB,sctrB) + B' * C * B * J1 * J2 * wt;
        end % end of loop on Gauss points
    end     % end of loop on elements  
end         % end of loop on patches

% Computing the external force vector
%%

% [W1,Q1] = quadrature(8, 'GAUSS', 1 );
% 
% % Loop over elements along left edge = noElemsV
% 
% elConnV   = patches(1).elConnV;
% noElemsV  = patches(1).noElemsV;
% 
% for e=1:noElemsV
%     xiE   = elRangeV(e,:); % [xi_i,xi_i+1]
%     conn  = elConnV(e,:);
%     pts   = leftPoints(conn,:);
%     sctrx = leftEdgeMesh(e,:);
%     sctry = sctrx + noCtrPts;
%     % loop over Gauss points
%     for gp=1:size(W1,1)
%         xi      = Q1(gp,:);
%         wt      = W1(gp);
%         Xi      = 0.5 * ( ( xiE(2) - xiE(1) ) * xi + xiE(2) + xiE(1));
%         J2      = 0.5 * ( xiE(2) - xiE(1) );
%         
%         [N dNdxi] = NURBS1DBasisDers(Xi,q,vKnot,weights);
%         
%         % compute the jacobian of physical and parameter domain mapping
%         % then the derivative w.r.t spatial physical coordinates
%         
%         x        = N     * pts; % global coord of GP
%         jacob1   = dNdxi * pts;
%         J1       = norm (jacob1);
%              
%         str      = exact_plate_hole(x,a);
%         tx       = -str(1);
%         ty       = -str(3);
%         fac      = J1 * J2 * wt;
%         f(sctrx) = f(sctrx) + N' * tx * fac;
%         f(sctry) = f(sctry) + N' * ty * fac;
%     end
% end
% 
% for e=1:noElemsV
%     xiE   = elRangeV(e,:); % [xi_i,xi_i+1]
%     conn  = elConnV(e,:);
%     pts   = rightPoints(conn,:);
%     sctrx = rightEdgeMesh(e,:);
%     sctry = sctrx + noCtrPts;
%     % loop over Gauss points
%     for gp=1:size(W1,1)
%         xi      = Q1(gp,:);
%         wt      = W1(gp);
%         Xi      = 0.5 * ( ( xiE(2) - xiE(1) ) * xi + xiE(2) + xiE(1));
%         J2      = 0.5 * ( xiE(2) - xiE(1) );
%         
%         [N dNdxi] = NURBS1DBasisDers(Xi,q,vKnot,weights);
%         
%         % compute the jacobian of physical and parameter domain mapping
%         % then the derivative w.r.t spatial physical coordinates
%         
%         x        = N     * pts; % global coord of GP
%         jacob1   = dNdxi * pts;
%         J1       = norm (jacob1);
%              
%         str      = exact_plate_hole(x,a);
%         tx       = str(1);
%         ty       = str(3);
%         fac      = J1 * J2 * wt;
%         f(sctrx) = f(sctrx) + N' * tx * fac;
%         f(sctry) = f(sctry) + N' * ty * fac;
%     end
% end

% for e=1:noElemsV
%     xiE   = elRangeV(e,:); % [xi_i,xi_i+1]
%     conn  = elConnV(e,:);
%     pts   = topPoints(conn,:);
%     sctrx = topEdgeMesh(e,:);
%     sctry = sctrx + noCtrPts;
%     % loop over Gauss points
%     for gp=1:size(W1,1)
%         xi      = Q1(gp,:);
%         wt      = W1(gp);
%         Xi      = 0.5 * ( ( xiE(2) - xiE(1) ) * xi + xiE(2) + xiE(1));
%         J2      = 0.5 * (   xiE(2) - xiE(1) );
%         
%         [N dNdxi] = NURBS1DBasisDers(Xi,q,vKnot,weights);
%         
%         % compute the jacobian of physical and parameter domain mapping
%         % then the derivative w.r.t spatial physical coordinates
%         
%         x        = N     * pts; % global coord of GP
%         jacob1   = dNdxi * pts;
%         J1       = norm (jacob1);
%              
%         str      = exact_plate_hole(x,a);
%         tx       = str(3);
%         ty       = str(2);
%         fac      = J1 * J2 * wt;
%         f(sctrx) = f(sctrx) + N' * tx * fac;
%         f(sctry) = f(sctry) + N' * ty * fac;
%     end
% end
% 
% for e=1:noElemsV
%     xiE   = elRangeV(e,:); % [xi_i,xi_i+1]
%     conn  = elConnV(e,:);
%     pts   = bottomPoints(conn,:);
%     sctrx = bottomEdgeMesh(e,:);
%     sctry = sctrx + noCtrPts;
%     % loop over Gauss points
%     for gp=1:size(W1,1)
%         xi      = Q1(gp,:);
%         wt      = W1(gp);
%         Xi      = 0.5 * ( ( xiE(2) - xiE(1) ) * xi + xiE(2) + xiE(1));
%         J2      = 0.5 * (   xiE(2) - xiE(1) );
%         
%         [N dNdxi] = NURBS1DBasisDers(Xi,q,vKnot,weights);
%         
%         % compute the jacobian of physical and parameter domain mapping
%         % then the derivative w.r.t spatial physical coordinates
%         
%         x        = N     * pts; % global coord of GP
%         jacob1   = dNdxi * pts;
%         J1       = norm (jacob1);
%              
%         str      = exact_plate_hole(x,a);
%         tx       = -str(3);
%         ty       = -str(2);
%         fac      = J1 * J2 * wt;
%         f(sctrx) = f(sctrx) + N' * tx * fac;
%         f(sctry) = f(sctry) + N' * ty * fac;
%     end
% end

disp([num2str(toc),'  APPLYING BOUNDARY  CONDITIONS'])

udofs=[leftNodes rightNodes]; % global indecies  of the fixed x displacements
vdofs=vNodes+noCtrPts; % global indecies  of the fixed y displacements

applyBC

% SOLVE SYSTEM

disp([num2str(toc),'  SOLVING THE SYSTEM'])
U=K\f;

Ux    = U(1:noCtrPts);
Uy    = U(noCtrPts+1:noDofs);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% POST-PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% Stress computation

vtuFile0 = 'bracket';

% Loop over patches

%figure
for ip=1:noPatches
    index      = patches(ip).index;
    elRangeU   = patches(ip).elRangeU;
    elRangeV   = patches(ip).elRangeV;
    element    = patches(ip).element;    
    elementL   = patches(ip).elementLocal;    
    uKnot      = patches(ip).uKnot;
    vKnot      = patches(ip).vKnot;
    controlPts = patches(ip).controlPts;
    weights    = patches(ip).weights;
    noPtsX     = patches(ip).noPtsX;
    noPtsY     = patches(ip).noPtsY;
    noElems    = patches(ip).noElemsU * patches(ip).noElemsV;
    
    vtuFile    = strcat(vtuFile0,num2str(ip));
    
    plotStressMP
end








