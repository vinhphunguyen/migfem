function [W,Q] = gaussianQuad2DNURBS(uOrder,vOrder)

% The function quadrature returns a n x 1 column vector W of quadrature


quadpoint  = zeros(uOrder*vOrder,2);
quadweight = zeros(uOrder*vOrder,1);

% Gauss quadrature along three directions

[ptU, wtU] = gaussQuad(uOrder);
[ptV, wtV] = gaussQuad(vOrder);

n=1;

for i = 1:uOrder
    for j = 1:vOrder
            quadpoint(n,:) = [ ptU(i), ptV(j)];
            quadweight(n)  = wtU(i)*wtV(j);
            n              = n+1;
    end
end

Q = quadpoint;
W = quadweight;

% END OF FUNCTION
