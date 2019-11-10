function [Xi]=global2LocalMap3D(gpoint,nodes)

% nodes
% gpoint

Xi    = zeros(1,3);
nMax  = 10;
tol   = 1e-14; tol=tol^2;

% shift coord. to center (better accuracy when mapping)
xm = sum(nodes,1)/size(nodes,1);

nodes(:,1) = nodes(:,1) - xm(1);
nodes(:,2) = nodes(:,2) - xm(2);
nodes(:,3) = nodes(:,3) - xm(3);

gpoint(:,1) =  gpoint(:,1) - xm(1);
gpoint(:,2) =  gpoint(:,2) - xm(2);
gpoint(:,3) =  gpoint(:,3) - xm(3);

% Newton-Raphson iterations
dSi = 1;
n   = 1;
while  dSi>tol && n<nMax
    [N,dNdxi]=lagrange_basis('B8',Xi);   % compute shape functions
    x       = N'*nodes;
    hessian = nodes' * dNdxi;
    f       = x - gpoint;
    dxi     = inv(hessian)*f';
    Xi      = Xi-dxi';
    n       = n+1;
    dSi     = dxi(1)^2+dxi(2)^2+dxi(3)^2;
%     dSi
%     f
end

if n==nMax && dSi>tol
    warning(['mapping Gauss points; residual, dX = ',num2str(sqrt(dSi))])
end