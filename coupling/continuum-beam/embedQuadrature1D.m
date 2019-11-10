function [W,Q] = embedQuadrature1D(  noGPs, elemState, x0, pts,e)

if elemState(e) == 1
    noGPs = noGPs*5;
end

[W,Q]=quadrature(  noGPs, 'GAUSS', 1 ); % 2 point quadrature
x=zeros(noGPs,1);

if elemState(e) == 1
    for gp=1:size(W,1)
        pt      = Q(gp,:);        
        [N,dNdxi]=lagrange_basis('L2',pt);  % element shape functions
        x= N'*pts(:,1);
    end
    id = find(x<x0);
    W(id) = [];
    Q(id) = [];
end