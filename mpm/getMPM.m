function [phi,dphi]=getMPM(x,h)

% compute the MPM (linear) shape function at point x
% h:  element size

if      ( abs(x) <= h )
    phi  = 1.0 - abs(x)/h;
    dphi = -1/h*sign(x);
else
    phi = 0;
    dphi = 0;
end