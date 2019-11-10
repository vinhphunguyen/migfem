function [K,fint] = getTangentStiffnessMatrix (U,C)
% compute the tangent stiffness matrix and the internal force
% vector for 2D Total Lagrangian formulation for finite deformation.
% VP Nguyen
% Cardiff University, March 2013.

global element p q index elRangeU elRangeV uKnot vKnot nElmLK...
    weights noCtrPts noDofs Ke0 controlPts W GP 

fint = zeros(noDofs,1);        % internal force vector
%K    = zeros(noDofs,noDofs);        % internal force vector

%% one the assembly
nElNod = size(element,2);
nElDof = nElNod*2;
nElmLK = nElDof^2;

jSprLK = 1:nElmLK;
vIdGrd = ones(1,nElDof);

noElems=size(element,1);
nSprGK = nElmLK*noElems;
vSprGK = zeros(nSprGK,1);
jSprRw = zeros(nSprGK,1);
jSprCl = zeros(nSprGK,1);

for e=1:noElems
    idu    = index(e,1);
    idv    = index(e,2);
    xiE    = elRangeU(idu,:); % [xi_i,xi_i+1]
    etaE   = elRangeV(idv,:); % [eta_j,eta_j+1]
    
    sctr   = element(e,:);          %  element scatter vector
    nn              = length(sctr);
    sctrB           = zeros(1,2*nn);
    sctrB(1:2:2*nn) = 2*sctr-1;
    sctrB(2:2:2*nn) = 2*sctr;
    pts    = controlPts(sctr,:);    %  element coordinates
    eDisp  = U(sctrB);              %  element displacement vector
    Ke     = Ke0;
    
    % loop over Gauss points
    
    for gp=1:size(W,1)
        pt      = GP(gp,:);
        wt      = W(gp);
        
        % compute coords in parameter space
        Xi      = parent2ParametricSpace(xiE,pt(1));
        Eta     = parent2ParametricSpace(etaE,pt(2));
        J2      = jacobianPaPaMapping(xiE,etaE);
        
        [dRdxi dRdeta] = NURBS2Dders([Xi; Eta],p,q,uKnot,vKnot,weights');
        
        % compute the jacobian of physical and parameter domain mapping
        % then the derivative w.r.t spatial physical coordinates
        
        jacob = [dRdxi; dRdeta]*pts;
        J1    = det(jacob);
        
        % Jacobian inverse and spatial derivatives
        
        invJacob   = inv(jacob);
        dRdx       = invJacob*[dRdxi;dRdeta];
        
        % compute gradient deformation and strain
        
        [strain,F] = getKinematics(dRdx, eDisp);
        
        % compute stress
        
        stress = C*strain;
        T      = stress2matrix( stress );
        
        % compute B matrices
        
        B     = getBmatrix  (dRdx, F);
        BNL   = getBNLmatrix(dRdx);
        
        % compute elementary stiffness matrix and internal force
        
        Ke          = Ke   + ( B' * C * B + BNL' * T * BNL ) * J1 * J2 * wt;
        fint(sctrB) = fint(sctrB) + B' * stress * J1 * J2 * wt;
    end
    
    %K(sctrB,sctrB) = K(sctrB,sctrB) + Ke;
    %----------------------------------------------------------------------
    % Global Assembly (Sparse)
    %----------------------------------------------------------------------
    
    vElDof = sctrB';
    mRwGrd = vElDof(:,vIdGrd);
    mClGrd = mRwGrd';
    
    jSprRw(jSprLK) = mRwGrd(:);
    jSprCl(jSprLK) = mClGrd(:);
    vSprGK(jSprLK) = Ke(:);
    
    jSprLK         = jSprLK + nElmLK; % move to the next position
end

% Here comes the total stiffness matrix in one shot!!!

K = sparse(jSprRw,jSprCl,vSprGK,noDofs,noDofs);

%% some utilities

function B=getBmatrix(dphi, F)
% compute the linear B matrix
% dphi: derivatives of shape functions wrt physical coords
% F:    gradient deformation matrix (2x2)

B = zeros(3, 2*length(dphi));

for i=1:length(dphi)
    B(1,2*i-1 )  = dphi(1,i)*F(1,1);
    B(1,2*i   )  = dphi(1,i)*F(2,1);
    B(2,2*i-1 )  = dphi(2,i)*F(1,2);
    B(2,2*i   )  = dphi(2,i)*F(2,2);
    B(3,2*i-1 )  = dphi(2,i)*F(1,1) + dphi(1,i)*F(1,2);
    B(3,2*i   )  = dphi(1,i)*F(2,2) + dphi(2,i)*F(2,1);
end

function Bnl = getBNLmatrix(dphi)
% compute the nonlinear B matrix
% dphi: derivatives of shape functions wrt physical coords

Bnl = zeros(4,2*length(dphi));

for i=1:length(dphi)
    Bnl(1,2*i-1) = dphi(1,i);
    Bnl(2,2*i-1) = dphi(2,i);
    Bnl(3,2*i)   = dphi(1,i);
    Bnl(4,2*i)   = dphi(2,i);
end

function T = stress2matrix( stress )
% Transform a vector stress to a matrix stress

T = zeros(4,4);

T(1,1) = stress(1);
T(2,2) = stress(2);
T(1,2) = stress(3);
T(2,1) = stress(3);

T(3:4,3:4) = T(1:2,1:2);

%%
function [strain,F] = getKinematics(dphi, elstate)
% compute the deformation gradient and the
% Green-Lagrange strain

% F_{ij} = I_{ij} + u_{i,j}
% F_{11} = 1+ u_1,1 = N_I,1 u_I1
% F_{12} =    u_1,2 = N_I,2 u_I1
% F_{21} =    u_2,1 = N_I,1 u_I2
% F_{22} = 1+ u_2,2 = N_I,2 u_I2

nn = length(dphi); % number of nodes per element

F = zeros(2);

F(1,1) = 1 + dphi(1,:) * elstate(1:2:2*nn);
F(1,2) =     dphi(2,:) * elstate(1:2:2*nn);
F(2,1) =     dphi(1,:) * elstate(2:2:2*nn);
F(2,2) = 1 + dphi(2,:) * elstate(2:2:2*nn);


E = 0.5*(F'*F-eye(2)); % E=1/2(F^T F - I)

strain(1,1) = E(1,1);
strain(2,1) = E(2,2);
strain(3,1) = 2.0*E(1,2);

%%

