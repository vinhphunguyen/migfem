function [f] = branchFunctions3d(Phi,N)
% Compute the branch functions spanning the near tip field for LEFM
% f   = f(r,theta)
% r   = r(phi,psi)
% phi = phi(xi,eta,zeta)
% Input:
%  Phi    : level sets at element under consideration
%  N      : shape functions
%  dNdxi  : derivatives of N w.r.t xi
% Vinh Phu Nguyen

% level sets at this point

normalPhi  = Phi(:,1);
tangentPhi = Phi(:,2);

phi        = N' * normalPhi;
psi        = N' * tangentPhi;

% (r,theta) coordinates

r     = sqrt(phi^2+psi^2);
theta = atan2(phi,psi);

if (r == 0) 
  disp('r=0')
end

r2   = sqrt(r);
fac  = 0.5/r2 ;

st2  = sin(theta/2.);
ct2  = cos(theta/2.);
st   = sin(theta)   ;
ct   = cos(theta)   ;

% ------------------
%     Functions 
% ------------------

f(1) = r2 * st2 ;
f(2) = r2 * ct2;
f(3) = r2 * st2 * st;
f(4) = r2 * ct2 * st;







