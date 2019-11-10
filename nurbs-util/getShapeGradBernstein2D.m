function [B dB] = getShapeGradBernstein2D(p,q,xi,eta)
%
% Compute Bernstein basis and first derivaties.
% The basis is defined in [-1,1] to mimic FEM.
%
% VP Nguyen, Cardiff University, UK

%% Initialization
ncpt  = (p+1)*(q+1);
B     = zeros(ncpt,1);
dB    = zeros(ncpt,2);

%% Calculation

[Bx dBdxi] = getShapeGradBernstein(p,xi);
[Be dBdet] = getShapeGradBernstein(q,eta);

for j=1:q+1
    for i=1:p+1
        id = (p+1)*(j-1)+i;
        B(id)    = Bx(i)    * Be(j);
        dB(id,1) = dBdxi(i) * Be(j);
        dB(id,2) = Bx(i)    * dBdet(j);
    end
end

