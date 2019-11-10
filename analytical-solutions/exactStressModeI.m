function [sigmax,sigmay,sigmaxy] = exactStressModeI(x,E,nu,sigmato,xTip,adv,cracklength)
% Compute the exact displacement for the infinite plate with centered crack
% Inputs:
%     - x(1,2) : coordinates of point where exact displacements are to be
%                evaluated
%     - E,nu: material properties
%     - sigmato : the loading
%     - xTip,adv,cracklength: crack data

alfa      = atan2(adv(2),adv(1)); %  inclination of local coord
QT        = [cos(alfa) sin(alfa); -sin(alfa) cos(alfa)];
xp        = QT*(x-xTip)';         % local coordinates
[theta,r] = cart2pol(xp(1),xp(2));% local polar coordinates

% if r==0
%     disp('Hang')
%     r = r +0.001;
% end

KI      = sigmato*sqrt(pi*cracklength);   % exact KI
KI      = 1;   % exact KI
facS    = KI/sqrt(2*pi*r);

sigmax  = facS*cos(theta/2)*(1-sin(theta/2)*sin(1.5*theta));
sigmay  = facS*cos(theta/2)*(1+sin(theta/2)*sin(1.5*theta));
sigmaxy = facS*cos(theta/2)*sin(theta/2)*cos(1.5*theta);





