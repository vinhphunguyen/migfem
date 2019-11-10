function [W,Q] = gaussianQuadNURBS(uOrder,vOrder,wOrder)

% The function quadrature returns a n x 1 column vector W of quadrature


quadpoint  = zeros(uOrder*vOrder*wOrder,3);
quadweight = zeros(uOrder*vOrder*wOrder,1);

% Gauss quadrature along three directions

[ptU, wtU] = gaussQuad(uOrder);
[ptV, wtV] = gaussQuad(vOrder);
[ptW, wtW] = gaussQuad(wOrder);

n=1;

for i = 1:uOrder
    for j = 1:vOrder
        for k = 1:wOrder
            quadpoint(n,:) = [ ptU(i), ptV(j), ptW(k) ];
            quadweight(n)  = wtU(i)*wtV(j)*wtW(k);
            n              = n+1;
        end
    end
end

Q = quadpoint;
W = quadweight;

% END OF FUNCTION
