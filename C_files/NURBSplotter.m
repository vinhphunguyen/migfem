% This routine uses the C-function to construct the B-spline functions
% Robert Simpson, Cardiff University, UK

close all
clear all
clc


knot=[0 0 0 0.25 0.5 0.75 1 1 1];
xi=0:0.01:max(knot);
numPts=numel(xi);
m = numel(knot)-1;
p = 2;
n = m - p;
weights = ones(1,n);
weights(2) = 1;

points = [0 0; 1 3; 2 2; 3 3; 4 5; 6 6; 7 8; 8 4];
interpolatedPoints = zeros(numPts,2);
weights(4) =1;
weights(6) = 1;

BsplineVals = zeros(numPts,n);
NURBSderivs = zeros(numPts,n);

tic
for i=1:n
    for c=1:numPts
        [BsplineVals(c,i+1) NURBSderivs(c,i+1)] = NURBSbasis(i, p, xi(c), knot, weights);
    end
end
toc

figure
plot(xi, BsplineVals)

figure
plot(xi, NURBSderivs)

for c=1:numPts
    [interpolatedPoints(c,1)] = NURBSinterpolation(xi(c), p, knot, points(:,1)', weights);
    [interpolatedPoints(c,2)] = NURBSinterpolation(xi(c), p, knot, points(:,2)', weights);      
end

figure
plot(interpolatedPoints(:,1), interpolatedPoints(:,2), 'k-', points(:,1), points(:,2), 'ko')
