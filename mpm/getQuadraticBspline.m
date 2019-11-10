function [phi,dphi]=getQuadraticBspline(x,h)

% compute the quadratic bspline at point x
% h:  element size
% derivatives: not yet implemented!!!


if      ( x < -1.5*h )
    phi  = 0;
    dphi = 0;
elseif      ( x >= -1.5*h ) && ( x <= -0.5*h )
    phi  = 0.5/h/h*x^2 + 1.5/h*x + 9/8;
    dphi = -1/h*sign(x);
elseif ( x <= 0.5*h)
    phi = -1/h/h*x^2 + 0.75;
    dphi = 0;
elseif ( x <= 1.5*h)
    phi = 0.5/h/h*x^2 -1.5/h*x + 9/8;
    dphi = 0;
else
    phi = 0;
    dphi = 0;
end