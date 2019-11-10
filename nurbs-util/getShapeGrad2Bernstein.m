function [B dBdxi dBdxi2]=getShapeGrad2Bernstein(p,xi)
%
% Compute Bernstein basis and first/second derivaties.
% The basis is defined in [-1,1] to mimic FEM.
%
% VP Nguyen, Cardiff University, UK

%% Initialization
ncpt   = p+1;
B      = zeros(ncpt,1);
dBdxi  = zeros(ncpt,1);
dBdxi2 = zeros(ncpt,1);

%% Calculation

for i=1:p+1
    B(i)      = bernstein(p,i,xi);
    dBdxi(i)  = 0.5*p*(bernstein(p-1,i-1,xi)-bernstein(p-1,i,xi));    
    dBdxi2(i) = 0.25*p*(p-1)*(bernstein(p-2,i-2,xi)-2*bernstein(p-2,i-1,xi)+bernstein(p-2,i,xi));  
end
