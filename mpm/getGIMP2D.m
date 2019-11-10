function [phi,dphi]=getGIMP2D(x,hx,hy,lpx,lpy)

% compute the 1D shape functions

[phiX,dphiX]=getGIMP(x(1),hx,lpx);
[phiY,dphiY]=getGIMP(x(2),hy,lpy);

% compute the 2D shape functions

phi     = phiX * phiY;
dphi(1) = dphiX* phiY;
dphi(2) =  phiX* dphiY;