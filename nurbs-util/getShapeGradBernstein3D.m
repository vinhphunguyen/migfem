function [B dB] = getShapeGradBernstein3D(p,q,r,xi,eta,zeta)
%
% Compute Bernstein basis and first derivaties.
% The basis is defined in [-1,1]x[-1,1]x[-1,1] to mimic FEM.
%
% VP Nguyen, Cardiff University, UK

%% Initialization
ncpt  = (p+1)*(q+1)*(r+1);
B     = zeros(ncpt,1);
dB    = zeros(ncpt,3);

%% Calculation

[Bx dBdxi] = getShapeGradBernstein(p,xi);
[Be dBdet] = getShapeGradBernstein(q,eta);
[Bz dBdez] = getShapeGradBernstein(r,zeta);

for k=1:r+1
    for j=1:q+1
        for i=1:p+1
            id = (p+1)*(j-1) + (p+1)*(q+1)*(k-1) + i;
            B(id)    = Bx(i)    * Be(j)    * Bz(k);
            dB(id,1) = dBdxi(i) * Be(j)    * Bz(k);
            dB(id,2) = Bx(i)    * dBdet(j) * Bz(k);
            dB(id,3) = Bx(i)    * Be(j)    * dBdez(k);
        end
    end
end

