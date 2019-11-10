function [ux,uy] = exactDispModeII(x,E,nu,stressState,sigmato,xTip,adv,cracklength)
% Compute the exact displacement for the infinite plate with centered crack
% under mode II loading
% Inputs:
%     - x(1,2) : coordinates of point where exact displacements are to be
%                evaluated
%     - E,nu, stressState: material properties
%     - sigmato : the loading
%     - xTip,adv,cracklength: crack data

alfa      = atan2(adv(2),adv(1)); %  inclination of local coord
QT        = [cos(alfa) sin(alfa); -sin(alfa) cos(alfa)];
xp        = QT*(x-xTip)';         % local coordinates
[theta,r] = cart2pol(xp(1),xp(2));% local polar coordinates

KII       = sigmato*sqrt(pi*cracklength);   % exact KII
KII       = 1;   % exact KII
mu        = E/(2*(1+nu));                   % shear modulus
fac       = 0.5*KII/mu*sqrt(r/2/pi);

if ( strcmp(stressState,'PLANE_STRAIN') )
    kappa = 3 - 4*nu;
else
    kappa = (3-nu)/(1+nu);
end

ux = fac*sin(theta/2)*( kappa + 1 + 2*cos(theta/2)^2);
uy = fac*cos(theta/2)*(-kappa + 1 + 2*sin(theta/2)^2);





