function [exact_stress] = exact_plate_hole(pt,a)
% Compute the exact stresses of the infinite plate with 
% a circular centered hole problem
% Vinh Phu Nguyen

x = pt(1,1);
y = pt(1,2);

r     = sqrt(x*x + y*y);
theta = atan(y/x);
%theta = atan2(y,x);

c2t = cos(2*theta);
c4t = cos(4*theta);
s2t = sin(2*theta);
s4t = sin(4*theta);
fac1 = (a/r)^2;
fac2 = fac1*fac1;
exact_stress(1) = 1-fac1*(1.5*c2t+c4t)+1.5*fac2*c4t;
exact_stress(2) =  -fac1*(0.5*c2t-c4t)-1.5*fac2*c4t;
exact_stress(3) =  -fac1*(0.5*s2t+s4t)+1.5*fac2*s4t;