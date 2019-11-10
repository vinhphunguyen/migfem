function [B dBdxi]=bernsteinBasis1D(p,xi)
%
% Compute univariate Bernstein basis and first derivaties.
% The basis is defined in [-1,1] to mimic FEM.
%
% VP Nguyen, Cardiff University, UK

%% Initialization
ncpt  = (p+1);
B     = zeros(ncpt,1);
dBdxi = zeros(ncpt,1);

%% Calculation

for i=1:p+1
    B(i)     = bernstein(p,i,xi);
    dBdxi(i) = 0.5*p*(bernstein(p-1,i-1,xi)-bernstein(p- 1,i,xi));
end


