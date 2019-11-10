function [sigmax,sigmay,sigmaxy] = exactStressModeII(x,E,nu,sigmato,xTip,adv,cracklength)
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

KII  = sigmato*sqrt(pi*cracklength);   % exact KII
KII  = 1;   % exact KII

facS = KII/sqrt(2*pi*r);
st2  = sin(theta/2);
ct2  = cos(theta/2);
c3t2 = cos(1.5*theta);

sigmax  = -facS*st2*(2 + ct2*c3t2);
sigmay  =  facS*st2*ct2*c3t2;
sigmaxy =  facS*ct2*(1-st2*sin(1.5*theta));





