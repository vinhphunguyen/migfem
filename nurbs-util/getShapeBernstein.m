function [B]=getShapeBernstein(p,xi)
%
% Compute Bernstein basis.
% The basis is defined in [-1,1] to mimic FEM.
%
% VP Nguyen, Cardiff University, UK

%% Initialization
ncpt  = p+1;
B     = zeros(ncpt,1);


%% Calculation

for i=1:p+1
    B(i)     = bernstein(p,i,xi);       
end
