function [f,dfdxi,dfdeta,dfdzeta] = branch3d(Phi,N,dNdxi,dNdeta,dNdzeta)
% Compute the branch functions spanning the near tip field for LEFM
% and its derivatives with respect to (xi,eta,zeta) coordinates.
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

% ---------------------------------------
%   Derivatives w.r.t (r,theta)
% ---------------------------------------

% first function

df1dr = fac * st2;
df1dt = 0.5 * r2 * ct2;

% second function

df2dr =  fac * ct2;
df2dt = -0.5 * fac * st2;

% third function

df3dr = fac * st2 * st;
df3dt = r2  * (0.5*ct2 * st + st2*ct);

% fourth function

df4dr = fac * ct2 * st;
df4dt = r2  * (-0.5*st2 * st + ct2*ct);

% -----------------------------------------
%   Derivatives of (r,theta) w.r.t phi, psi
% -----------------------------------------

drdphi = phi/r;
drdpsi = psi/r;
dtdphi = psi/r/r;
dtdpsi = -phi/r/r;

% ----------------------------------------------
%   Derivatives of (phi, psi) w.r.t xi,eta,zeta
% ----------------------------------------------

dphidxi   = dNdxi'   * normalPhi;
dphideta  = dNdeta'  * normalPhi;
dphidzeta = dNdzeta' * normalPhi;

dpsidxi   = dNdxi'   * tangentPhi;
dpsideta  = dNdeta'  * tangentPhi;
dpsidzeta = dNdzeta' * tangentPhi;

% ----------------------------------------------
%  Finally, derivatives of branch functions 
%  w.r.t (xi,eta,zeta)
% ----------------------------------------------

drdxi   = drdphi * dphidxi   + drdpsi * dpsidxi;
drdeta  = drdphi * dphideta  + drdpsi * dpsideta;
drdzeta = drdphi * dphidzeta + drdpsi * dpsidzeta;

dtdxi   = dtdphi * dphidxi   + dtdpsi * dpsidxi;
dtdeta  = dtdphi * dphideta  + dtdpsi * dpsideta;
dtdzeta = dtdphi * dphidzeta + dtdpsi * dpsidzeta;

% xi direction

dfdxi(1)   = df1dr * drdxi + df1dt * dtdxi;
dfdxi(2)   = df2dr * drdxi + df2dt * dtdxi;
dfdxi(3)   = df3dr * drdxi + df3dt * dtdxi;
dfdxi(4)   = df4dr * drdxi + df4dt * dtdxi;

% eta direction

dfdeta(1)  = df1dr * drdeta + df1dt * dtdeta;
dfdeta(2)  = df2dr * drdeta + df2dt * dtdeta;
dfdeta(3)  = df3dr * drdeta + df3dt * dtdeta;
dfdeta(4)  = df4dr * drdeta + df4dt * dtdeta;

% zeta direction

dfdzeta(1) = df1dr * drdzeta + df1dt * dtdzeta;
dfdzeta(2) = df2dr * drdzeta + df2dt * dtdzeta;
dfdzeta(3) = df3dr * drdzeta + df3dt * dtdzeta;
dfdzeta(4) = df4dr * drdzeta + df4dt * dtdzeta;





