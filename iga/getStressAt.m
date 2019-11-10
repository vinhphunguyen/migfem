function [sigma] = getStressAt(e,u)

global index elRangeU elRangeV element controlPts W Q C

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

pt      = Q(1,:);
wt      = W(1);

[R,dRdx,J1] = getShapeGrads2D(pt,xiE,etaE,pts);
B           = getBmatrix2D(B,dRdx); % B matrix
sigma       = C*B*u(sctrB);