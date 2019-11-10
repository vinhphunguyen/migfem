function [ux,uy,uz] = exactGriffith3D(phi,psi,E,nu,sigmato,cracklength)
% Compute the exact displacement for the infinite plate with centered crack
% Inputs:
%     - x(1,2) : coordinates of point where exact displacements are to be
%                evaluated
%     - E,nu: material properties
%     - sigmato : the loading
%     - xTip,adv,cracklength: crack data

r     = sqrt(phi^2+psi^2);
theta = atan2(phi,psi);

sigma = sigmato;
KI    = sigmato*sqrt(pi*cracklength);   % exact KI


mu=E/(2*(1+nu));
fac  = 2*(1+nu)*KI/E*sqrt(r/2/pi);

ux = fac*cos(theta/2)*(2-2*nu-cos(theta/2)^2);
uz = fac*sin(theta/2)*(2-2*nu-cos(theta/2)^2);
uy = 0.;





