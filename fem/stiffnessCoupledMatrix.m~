function [Ke,Ce,sctrB,sctrLB] = stiffnessCoupledMatrix...
                                (e,element1,element2,node1,node2,W,Q)

global elemType numdofs C alpha1 alpha2 ne lagDom l couDom2 fine2VirtualMap
global E0

sctr   = element1(e,:);
sctrB  = [ sctr sctr+numdofs ];

% Old code

% if e <= 40
%     sctrL  = lagDom(1,:);
%     sctrLL = element2(couDom2(1),:);
%     sctrLB = [ sctrL sctrL+numdofs ];
% else
%     sctrL  = lagDom(2,:);
%     sctrLL = element2(couDom2(2),:);
%     sctrLB = [ sctrL sctrL+numdofs ];
% end

i = fine2VirtualMap(e);

sctrL  = lagDom(i,:);
sctrLL = element2(couDom2(i),:);
sctrLB = [ sctrL sctrL+numdofs ];


nodesV = node2(sctrLL,:);

nn=length(sctr);

Ke = zeros(8,8);
Ce = zeros(8,8);

for q=1:size(W,1)
    pt=Q(q,:);
    wt=W(q);
    [N,dNdxi]=lagrange_basis(elemType,pt); % element shape functions
    J00    = node1(sctr,:)'*dNdxi;          % element Jacobian matrix
    invJ0 = inv(J00);
    dNdx  = dNdxi*invJ0;
    
    % for fine elements
    
    B(1,1:nn)       = dNdx(:,1)';
    B(2,nn+1:2*nn)  = dNdx(:,2)';
    B(3,1:nn)       = dNdx(:,2)';
    B(3,nn+1:2*nn)  = dNdx(:,1)';
    
    Nm(1,1:nn)       = N';
    Nm(2,nn+1:2*nn)  = N';
    
    % for coase element
    
    gpoint = N'*node1(sctr,:);
    
    lpoint = global2LocalQ4(nodesV,gpoint);
    
    [N,dNdxi]= lagrange_basis(elemType,lpoint);
    J0       = nodesV'*dNdxi;
    invJ0    = inv(J0);
    dNdx     = dNdxi*invJ0;
    
    BV(1,1:nn)       = dNdx(:,1)';
    BV(2,nn+1:2*nn)  = dNdx(:,2)';
    BV(3,1:nn)       = dNdx(:,2)';
    BV(3,nn+1:2*nn)  = dNdx(:,1)';
    
    NV(1,1:nn)       = N';
    NV(2,nn+1:2*nn)  = N';
    
    % COMPUTE ELEMENT STIFFNESS AT QUADRATURE POINT
    
    Ke = Ke + alpha1*BV'*C*B*W(q)*det(J00);
    Ce = Ce + NV'*Nm*W(q)*det(J00) + BV'*B*W(q)*det(J00)*l*l;
    %Ce = Ce + E0*NV'*Nm*W(q)*det(J00);
end