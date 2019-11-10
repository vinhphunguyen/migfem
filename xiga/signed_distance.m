% **********************************************************
%       TWO DIMENSIONAL ELEMENT FREE GALERKIN CODE
%                 Nguyen Vinh Phu
%            LTDS, ENISE, Juillet 2006
% *********************************************************

function d = signed_distance(xCr,pt)
% Compute the signed distance from point x to crack xCr
% Inputs:
%     - xcR(2,2) : coordinates of points defining the crack
%     - pt(1,2)   : coordinate of point 


x0  = xCr(1,1); y0 = xCr(1,2);
x1  = xCr(2,1); y1 = xCr(2,2);

x   = pt(1,1);
y   = pt(1,2);
l   = sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)) ;
phi = (y0-y1)*x + (x1-x0)*y + (x0*y1-x1*y0) ;
d   = phi/l;            

