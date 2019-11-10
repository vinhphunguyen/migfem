function [BsplineVals, NURBSderivs,xi] = getAllShapeGrads (uKnot,p,weights)

xi     = 0:0.01:max(uKnot);
numPts = numel(xi);
m      = numel(uKnot)-1;
n      = m - p;

BsplineVals = zeros(numPts,n);
NURBSderivs = zeros(numPts,n);

for i=1:n
    for c=1:numPts
        [BsplineVals(c,i+1) NURBSderivs(c,i+1)] = ...
            NURBSbasis(i, p, xi(c), uKnot, weights);
    end
end