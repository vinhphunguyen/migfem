function [Ke,sctrB] = stiffnessMatrix (e,element1,node1,W,Q)

global elemType numdofs C

sctr  = element1(e,:);           % element scatter vector
sctrB = [sctr sctr+numdofs]; % vector that scatters a B matrix
nn    = length(sctr);

Ke = zeros(8,8);

for q=1:size(W,1)
    pt = Q(q,:);
    wt = W(q);
    [N,dNdxi]=lagrange_basis(elemType,pt); % element shape functions
    J0    = node1(sctr,:)'*dNdxi;                % element Jacobian matrix
    invJ0 = inv(J0);
    dNdx  = dNdxi*invJ0;
    
    B(1,1:nn)       = dNdx(:,1)';
    B(2,nn+1:2*nn)  = dNdx(:,2)';
    B(3,1:nn)       = dNdx(:,2)';
    B(3,nn+1:2*nn)  = dNdx(:,1)';
    
    % COMPUTE ELEMENT STIFFNESS AT QUADRATURE POINT
    Ke=Ke+B'*C*B*W(q)*det(J0);
end  % of quadrature loop