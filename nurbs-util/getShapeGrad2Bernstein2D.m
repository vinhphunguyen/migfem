function [B dB dB2]=getShapeGrad2Bernstein2D(p,q,xi,eta)
%
% Compute Bernstein basis and first and second derivaties.
% The basis is defined in [-1,1] to mimic FEM.
%
% VP Nguyen, Cardiff University, UK

%% Initialization
ncpt   = (p+1)*(q+1);
B      = zeros(ncpt,1);
dB     = zeros(ncpt,2);
dB2    = zeros(ncpt,3);

%% Calculation

[Bx dBdxi dBdxi2] = getShapeGrad2Bernstein(p,xi);
[Be dBdet dBdet2] = getShapeGrad2Bernstein(q,eta);

for j=1:q+1
    for i=1:p+1
        id        = (p+1)*(j-1)+i;
        B(id)     = Bx(i)    * Be(j);
        
        dB(id,1)  = dBdxi(i) * Be(j);
        dB(id,2)  = Bx(i)    * dBdet(j);
        
        dB2(id,1) = dBdxi2(i) * Be(j);
        dB2(id,2) = dBdet2(j) * Bx(i);
        dB2(id,3) = dBdxi (i) * dBdet(j);
    end
end

