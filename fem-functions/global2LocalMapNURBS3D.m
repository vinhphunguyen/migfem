function [Xi]=global2LocalMapNURBS3D(x,xn,xiE,etaE,zetaE,p,q,r,uKnot,vKnot,wKnot,weights)
%
% Inverse mapping for 3D NURBS
% Given global point denoted by x, compute Xi which is the parameter space
% that corresponds to x.
% xn:              control points define element under consideration
% xiE,etaE,zetaE:  paramter ranges of element under consideration
% other params: evident
% Written by
% Vinh Phu Nguyen
% Cardiff University, Wales, UK
% June 2013

% parameters relating to the Newton-Raphson iterative procedure
nMax  = 10;
tol   = 1e-12; tol=tol^2;

% shift coord. to center (better accuracy when mapping)
% xm = sum(xn,1)/size(xn,1);
% 
% xn(:,1) = xn(:,1) - xm(1);
% xn(:,2) = xn(:,2) - xm(2);
% xn(:,3) = xn(:,3) - xm(3);
% 
% x(:,1) =  x(:,1) - xm(1);
% x(:,2) =  x(:,2) - xm(2);
% x(:,3) =  x(:,3) - xm(3);

% initial guess is the center of the element 
Xi_zrs = [sum(xiE)/2 sum(etaE)/2 sum(zetaE)/2];
Xi = Xi_zrs; 

% Newton-Raphson iterations
dSi = 1;
n   = 1;
while  dSi>tol && n<nMax
    [N dRdxi dRdeta dRdzeta] = NURBS3DBasisDers(Xi,p,q,r,uKnot,vKnot,wKnot,weights');
    xk      = N*xn;
    dNdxi   = [dRdxi' dRdeta' dRdzeta'];
    hessian = xn' * dNdxi;
    f       = xk - x;
    dxi     = inv(hessian)*f';
    Xi      = Xi-dxi';
    n       = n+1;
    dSi     = dxi(1)^2+dxi(2)^2+dxi(3)^2;

end

%    if Xi(1) < 0 || Xi(1) > 1  || ...
%       Xi(2) < 0 || Xi(2) > 1  || ...
%       Xi(3) < 0 || Xi(3) > 1  
%       error('error')
%    end
   
if n==nMax && dSi>tol
    warning(['mapping Gauss points; residual, dX = ',num2str(sqrt(dSi))])
end