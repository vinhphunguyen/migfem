function [f] = branchNode(phi,psi)
% Compute the branch functions spanning the near tip field for LEFM
% Inputs: 
% phi: normal level set 
% psi: tangent level set

r     = sqrt(phi^2+psi^2);
theta = atan2(phi,psi);

if( r ~=0 )
    r2 = sqrt(r);
else
    r2    = 0.1d-4;
    theta = 0.0d0 ;
end

st2  = sin(theta/2.);
ct2  = cos(theta/2.);
st   = sin(theta);
ct   = cos(theta);

% Functions 

f(1) = r2 * st2 ;
f(2) = r2 * ct2;
f(3) = r2 * st2 * st;
f(4) = r2 * ct2 * st;


        
