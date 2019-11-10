function [K,M,Cc] = getMatricesWithLumping ()

%
% Compute stiffness matrix and mass matrix for 2D continuum elements.
%
% VP Nguyen
% Cardiff University

global element controlPts p q uKnot vKnot weights index elRangeU elRangeV
global Q W rho C noDofs noCtrPts Ke0 Me0 damping

noElems = size(element,1);

%% one the assembly
nElNod = size(element,2);
nElDof = nElNod*2;
nElmLK = nElDof^2;
nSprGK = nElmLK*noElems;

jSprLK = 1:nElmLK;
vIdGrd = ones(1,nElDof);

% values, row indices, columns indices of the global K matrix
vSprGK1 = zeros(nSprGK,1);
vSprGK2 = zeros(nSprGK,1);
jSprRw  = zeros(nSprGK,1);
jSprCl  = zeros(nSprGK,1);

Ke = Ke0;
Me = Ke0;
Mm = zeros(nElNod,nElNod);

K    = zeros(noDofs,noDofs); % element Ke
M    = zeros(noDofs,1); % element Ke

% Loop over elements (knot spans)

for e=1:noElems
    idu    = index(e,1);
    idv    = index(e,2);
    xiE    = elRangeU(idu,:); % [xi_i,xi_i+1]
    etaE   = elRangeV(idv,:); % [eta_j,eta_j+1]
    
    sctr   = element(e,:);          %  element scatter vector
    nn     = length(sctr);
    
    sctrB(1,1:2:2*nn) = 2*sctr-1;    
    sctrB(1,2:2:2*nn) = 2*sctr  ;
    
    pts    = controlPts(sctr,:);
    
    B      = zeros(3,2*nn);
    N      = zeros(2,2*nn);
    
    % loop over Gauss points
    
    for gp=1:size(W,1)
        pt      = Q(gp,:);
        wt      = W(gp);
        
        [R,dRdx,J]  = getShapeGrads2D(pt,xiE,etaE,pts);                       
        B           = getBmatrix2D(B,dRdx); 
                
        % compute elementary stiffness matrix and
        % assemble it to the global matrix
        
        K(sctrB,sctrB) = K(sctrB,sctrB) + B' * C * B * J * wt;
        
        mQPt= rho*R'*R*J*wt;
        mQPt=sum(mQPt)';
        M(2*sctr-1)           = M(2*sctr-1)+mQPt;
        M(2*sctr)   = M(2*sctr)+mQPt;        
    end

    
    vElDof = sctrB';
    mRwGrd = vElDof(:,vIdGrd);
    mClGrd = mRwGrd';
    
    jSprRw(jSprLK)  = mRwGrd(:);
    jSprCl(jSprLK)  = mClGrd(:);
    vSprGK1(jSprLK) = Ke(:);
    vSprGK2(jSprLK) = Me(:);
    jSprLK          = jSprLK + nElmLK; % move to the next position
end

% Here comes the total stiffness matrix in one shot!!!

% K = sparse(jSprRw,jSprCl,vSprGK1,noDofs,noDofs);
% M = sparse(jSprRw,jSprCl,vSprGK2,noDofs,noDofs);
Cc = M * damping/rho;