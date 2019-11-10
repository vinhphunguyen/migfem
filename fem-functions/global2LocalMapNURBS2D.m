function X = global2LocalMapNURBS2D(x,xn,xiE,etaE,p, q, uKnot, vKnot, weights)
%==========================================================================
% Written by: Danas Stuala, Cardiff University, Wales, UK
% Adapted by: Vinh Phu Nguyen, Cardiff University, Wales, UK
%
% Inputs:
%  x : points in global coordinate system nx x 2 matrix
% xn : nodal coordinates of the element; nn x 2 matrix
%      where nn denotes the number of nodes.
% elemType: element type, 'Q4' or 'T3' etc.

nx = size(x,1);
nn = size(xn,1);

if ~(nx && nn)
    X = []; return
end

% control
n_max = 10; tol = 1e-14; %tol = tol^2; % dXi <= 1e-12

% shift coord. to center (better accuracy when mapping)
% xm = sum(xn,1)/nn;
% 
% xn(:,1) = xn(:,1) - xm(1);
% xn(:,2) = xn(:,2) - xm(2);
% 
% x(:,1) =  x(:,1) - xm(1);
% x(:,2) =  x(:,2) - xm(2);

X(nx,2) = 0; dXdx(2,2) = 0; Xi_zrs(1,2) = 0;

Xi_zrs = [sum(xiE)/2 sum(etaE)/2];

for ix = 1:nx
    
    n = 0; xi = x(ix,:); Xi = Xi_zrs; dSi = 1;
    
    while  dSi>tol && n<n_max                
%         Xii      = parent2ParametricSpace(xiE,Xi(1));
%         Eta     = parent2ParametricSpace(etaE,Xi(2));
                
        [R dRdxi dRdeta] = NURBS2DBasisDers(Xi,p,q,uKnot,vKnot,weights');
        dNdX = [dRdxi; dRdeta];
        dxdX = dNdX*xn;
        detJ = dxdX(1)*dxdX(4)-dxdX(2)*dxdX(3);
        
        dXdx(1) =  dxdX(4);
        dXdx(2) = -dxdX(2);
        dXdx(3) = -dxdX(3);
        dXdx(4) =  dxdX(1);
        
        dXi = detJ\(xi-R*xn)*dXdx;
        dSi = dXi(1)^2+dXi(2)^2;
        
        Xi = Xi+dXi; n = n+1;        
    end
    
    if n==n_max && dSi>tol
        warning(['mapping Gauss points; residual, dX = ',num2str(sqrt(dSi))])
    end
    
    X(ix,:) = Xi;    
end

%==========================================================================
end
