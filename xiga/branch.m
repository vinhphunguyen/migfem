% *************************************************************************
%                 TWO DIMENSIONAL ELEMENT FREE GALERKIN CODE
%                            Nguyen Vinh Phu
%                        LTDS, ENISE, Juillet 2006
% *************************************************************************
function [f,dfdx,dfdy] = branch(r,theta,alpha)
% Compute the branch functions spanning the near tip field for LEFM
% Inputs: 
%   (r,theta) : polar coordinates of points where the branch
%               functions are to be evaluated
%   alpha     : inclination of the crack tip segment w.r.t x axis

if( r ~=0 )
    r2 = sqrt(r);
else
    r2    = 0.1d-4;
    theta = 0.0d0 ;
end

fac  = 0.5/r2 ;

st2  = sin(theta/2.)    ;
ct2  = cos(theta/2.)    ;
s3t2 = sin(1.5 * theta) ;
c3t2 = cos(1.5 * theta) ;
st   = sin(theta)       ;
ct   = cos(theta)       ;

% ------------------
%     Functions 
% ------------------

f(1) = r2 * st2 ;
f(2) = r2 * ct2;
f(3) = r2 * st2 * st;
f(4) = r2 * ct2 * st;

% ---------------------------------------
%   Derivatives in local coords (x1,x2)
% ---------------------------------------

% first function

dPhi1dx1 = -fac * st2;
dPhi1dx2 = fac * ct2;


% second function

dPhi2dx1  = dPhi1dx2;
dPhi2dx2  = -dPhi1dx1;


% third function

dPhi3dx1 = -fac * s3t2 * st ;
dPhi3dx2 = fac * (st2 + s3t2 * ct);


% fourth function

dPhi4dx1  = -fac * c3t2 * st ;
dPhi4dx2  =  fac * (ct2 + c3t2 * ct);


% ---------------------------------------
%   Derivatives in global coords (x,y)
% ---------------------------------------

dx1dx = cos(alpha); dx2dx = -sin(alpha);
dx1dy = sin(alpha); dx2dy = cos(alpha);

dfdx(1) = dPhi1dx1 * dx1dx + dPhi1dx2 * dx2dx  ;
dfdy(1) = dPhi1dx1 * dx1dy + dPhi1dx2 * dx2dy  ;

dfdx(2) = dPhi2dx1 * dx1dx + dPhi2dx2 * dx2dx  ;
dfdy(2) = dPhi2dx1 * dx1dy + dPhi2dx2 * dx2dy  ;

dfdx(3) = dPhi3dx1 * dx1dx + dPhi3dx2 * dx2dx  ;
dfdy(3) = dPhi3dx1 * dx1dy + dPhi3dx2 * dx2dy  ;

dfdx(4) = dPhi4dx1 * dx1dx + dPhi4dx2 * dx2dx  ;
dfdy(4) = dPhi4dx1 * dx1dy + dPhi4dx2 * dx2dy  ;
