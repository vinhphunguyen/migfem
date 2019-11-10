% This script plot H(xi)*R(xi) function where
% H(xi) is the Heaviside function and R(xi) are
% the NURBS basis

addpath('~/code/xfem-efg-matlab/fem_util');
addpath C_files/

close all
clear all
clc


knot=[0 0 0 0.25 0.5 0.75 1 1 1];
xi=0:0.005:max(knot);
numPts=numel(xi);
m = numel(knot)-1;
p = 2;
n = m - p;
weights = ones(1,n);

% compute NURBS basis 

BsplineVals1 = zeros(numPts,n);
BsplineVals2 = zeros(numPts,n);

NURBSderivs = zeros(numPts,n);

for i=1:n
    for c=1:numPts
        [BsplineVals1(c,i) NURBSderivs(c,i)] = ...
            NURBSbasis(i, p, xi(c), knot, weights);
    end
end

% plot NURBS basis

figure
plot(xi, BsplineVals1,'LineWidth',1.2)

% compute discontinuous NURBS basis

xiBar = 0.4;
uspan = FindSpan(n-1,p,xiBar,knot);


for i=1:p+1
  cutId(i) = uspan-2+i;
end

for i=1:n
    for c=1:numPts
        dist = -xi(c)+xiBar;
        H    = heaviside(dist);
        if find(cutId==i)
            BsplineVals2(c,i)=BsplineVals1(c,i)*H;           
        else
            BsplineVals2(c,i)=BsplineVals1(c,i);           
        end
    end
end

figure
plot(xi, BsplineVals2,'LineWidth',1.2)
