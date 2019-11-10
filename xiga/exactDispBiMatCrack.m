function [ux,uy] = exactDispBiMatCrack(r,theta,e,g1,k1)

% Compute the near tip displacement field of an interface crack
% (r,theta): polar coordinates
% e,g1,k1: Dundurs parameter, shear modulus and Kosolov constant of mat1
% Vinh Phu Nguyen
% Hue, Vietnam April 2012

% upper half plane (e*pi) and lower half plane (-e*pi)

if ((theta >= 0) && (theta < pi)) ||...
        ((theta < -pi) && (theta >= -2*pi))
    epi = e*pi;                                   % Constant, delta, upper half plane
elseif ((theta >= pi) && (theta < 2*pi)) ||...
        ((theta < 0) && (theta >= -pi))
    epi = -e*pi;                                    % Constant, delta, lower half plane
end

A = exp(-epi+e*theta)/(1+4*e^2)/cosh(e*pi);

st  = sin(theta);
ct  = cos(theta);
ct2 = cos(theta/2);
st2 = sin(theta/2);

a1  = cos(e*log(r));
a2  = sin(e*log(r));

fac = 1+4*e^2;

K1  = 1;
K2  = 1;

u1I  = A*(-exp(2*epi-2*e*theta)*(ct2 + 2*e*st2) + k1*(ct2 - 2*e*st2) + fac*st2*st);
u1II = A*( exp(2*epi-2*e*theta)*(st2 - 2*e*ct2) + k1*(st2 + 2*e*ct2) + fac*ct2*st);
u2I  = A*( exp(2*epi-2*e*theta)*(st2 - 2*e*ct2) + k1*(st2 + 2*e*ct2) - fac*ct2*st);
u2II = A*( exp(2*epi-2*e*theta)*(ct2 + 2*e*st2) - k1*(ct2 - 2*e*st2) + fac*st2*st);

fac  = (0.5/g1)*sqrt(r/2/pi);

ux = fac*( (K1*a1 - K2*a2)*u1I + (K1*a2 + K2*a1)*u1II );
uy = fac*( (K1*a1 - K2*a2)*u2I + (K1*a2 + K2*a1)*u2II );










