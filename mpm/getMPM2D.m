function [phi,dphi]=getMPM2D(x,hx,hy)

% compute the 1D shape functions

[phiX,dphiX]=getMPM(x(1),hx);
[phiY,dphiY]=getMPM(x(2),hy);

% compute the 2D shape functions

phi     = phiX * phiY;
dphi(1) = dphiX* phiY;
dphi(2) =  phiX* dphiY;