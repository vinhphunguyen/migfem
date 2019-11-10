function [exact_stress disp] = exact_solution_hole(pt,a,E,nu)
% Compute the exact stresses of the infinite plate with 
% a circular centered hole problem

x     = pt(1,1);
y     = pt(1,2);

r     = sqrt(x*x + y*y);
theta = atan(y/x);

c1t = cos(theta);
c2t = cos(2*theta);
c3t = cos(3*theta);
c4t = cos(4*theta);
s1t = sin(theta);
s2t = sin(2*theta);
s3t = sin(3*theta);
s4t = sin(4*theta);
fac1 = (a/r)^2;
fac2 = fac1*fac1;

exact_stress(1) = 1-fac1*(1.5*c2t+c4t)+1.5*fac2*c4t;
exact_stress(2) =  -fac1*(0.5*c2t-c4t)-1.5*fac2*c4t;
exact_stress(3) =  -fac1*(0.5*s2t+s4t)+1.5*fac2*s4t;

mu   = E/2/(1+nu); % shear modulus
fac  = a/(8*mu);
fac1 = 2*(a/r)^3;
fac2 = 2*a/r;

%k    = 3-4*nu; % plane strain
k    = (3-nu)/(1+nu);

disp(1) = -fac*( r/a*(k+1)*c1t + fac2*( (k+1)*c1t + c3t ) - fac1*c3t );
disp(2) = -fac*( r/a*(k-3)*s1t + fac2*( (1-k)*s1t + s3t ) - fac1*s3t );



