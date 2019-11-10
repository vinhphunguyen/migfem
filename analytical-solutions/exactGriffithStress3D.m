
function [sigmax,sigmay,sigmaz,sigmayz,sigmaxz,sigmaxy] = exactGriffithStress3D...
    (phi,psi,E,nu,sigmato,cracklength)
% Compute the exact displacement for the infinite plate with centered crack
% Inputs:
%     - phi, psi: level set values
%     - E,nu: material properties
%     - sigmato : the loading
%     - cracklength: crack data

r     = sqrt(phi^2+psi^2);
theta = atan2(phi,psi);

sigma = sigmato;
KI    = sigmato*sqrt(pi*cracklength);   % exact KI

mu   =E/(2*(1+nu));
facS = KI/sqrt(2*pi*r);


sigmax  = facS*cos(theta/2)*(1-sin(theta/2)*sin(1.5*theta));
sigmaz  = facS*cos(theta/2)*(1+sin(theta/2)*sin(1.5*theta));
sigmaxz = facS*cos(theta/2)*sin(theta/2)*cos(1.5*theta);
sigmay  = nu*(sigmax+sigmaz);
sigmaxy = 0;
sigmayz = 0;




