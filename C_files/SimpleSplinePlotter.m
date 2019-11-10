% Robert Simpson 
% Cardiff University, UK

close all
clear all
clc


knot=[0 0 0 1 2 3 4 4 4]; 
xi=0:0.1:max(knot);
numPts=numel(xi);
m = numel(knot)-1;
p = 2;
n = m - p;
weights = ones(1,n);
weights(2) = 1;

points = [0 0; 0.5 0.5; 1 1.5; 2 2; 3 0.5; 4 2];
interpolatedPoints = zeros(numPts,2);
% weights(2) =2;
weights(2) =1.5;
weights(3) =2.7;

BsplineVals = zeros(numPts,n);

tic
for i=1:n
    for c=1:numPts
        [BsplineVals(c,i+1) ~] = NURBSbasis(i, p, xi(c), knot, weights);
    end
end
toc

figure(1)
plot(xi, BsplineVals)

%% and now plot the spline

for c=1:numPts
    [interpolatedPoints(c,1)] = NURBSinterpolation(xi(c), p, knot, points(:,1)', weights);
    [interpolatedPoints(c,2)] = NURBSinterpolation(xi(c), p, knot, points(:,2)', weights);
end

figure(3); hold on
plot(interpolatedPoints(:,1), interpolatedPoints(:,2), 'k-', points(:,1), points(:,2), 'ko')
hold off

%% and save data for plot program

save 'dat_files/coords' points -ASCII
save 'dat_files/interpolatedPoints.dat' interpolatedPoints -ASCII

