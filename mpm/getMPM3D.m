function [phi,dphi]=getMPM3D(x,hx,hy,hz)

% compute the 1D shape functions

[phiX,dphiX]=getMPM(x(1),hx);
[phiY,dphiY]=getMPM(x(2),hy);
[phiZ,dphiZ]=getMPM(x(3),hz);

% compute the 2D shape functions

phi     = phiX  * phiY * phiZ;
dphi(1) = dphiX * phiY * phiZ;
dphi(2) =  phiX * dphiY* phiZ;
dphi(3) =  phiX * phiY * dphiZ;