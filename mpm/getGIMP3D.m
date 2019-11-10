function [phi,dphi]=getGIMP3D(x,hx,hy,hz,lpx,lpy,lpz)

% compute the 1D shape functions

[phiX,dphiX]=getGIMP(x(1),hx,lpx);
[phiY,dphiY]=getGIMP(x(2),hy,lpy);
[phiZ,dphiZ]=getGIMP(x(3),hz,lpz);

% compute the 2D shape functions

phi     = phiX  * phiY * phiZ;
dphi(1) = dphiX * phiY * phiZ;
dphi(2) =  phiX * dphiY* phiZ;
dphi(3) =  phiX * phiY * dphiZ;