addpath nurbs/
addpath C_files/
addpath('~/code/xfem-efg-matlab/fem_util');

% replicate the example in Figure 2.23 of the IGA book
% Cottrel, Hughes and Bazilevs

%% original curve (one linear element)

controlPts = zeros(2,3);

controlPts = [0 1;
              0 0];
          
uKnot       = [0 0 1 1];          
linearCurve = nrbmak(controlPts, uKnot);


% plot basis

xi=0:0.01:max(uKnot);
numPts=numel(xi);
m = numel(uKnot)-1;
p = linearCurve.order-1;
n = m - p;
weights = ones(1,n);

BsplineVals = zeros(numPts,n);
NURBSderivs = zeros(numPts,n);

for i=1:n
    for c=1:numPts
        [BsplineVals(c,i+1) NURBSderivs(c,i+1)] = ...
            NURBSbasis(i, p, xi(c), uKnot, weights);
    end
end

figure
plot(xi, BsplineVals,'LineWidth',1.4)

%% knot insertion

newKnots = [1/3 2/3];

refineLinearCurve = nrbkntins(linearCurve,newKnots);
% plot basis

uKnot = refineLinearCurve.knots;
xi=0:0.001:max(uKnot);
numPts=numel(xi);
m = numel(uKnot)-1;
p = refineLinearCurve.order-1;
n = m - p;
weights = ones(1,n);

BsplineVals = zeros(numPts,n);
NURBSderivs = zeros(numPts,n);

for i=1:n
    for c=1:numPts
        [BsplineVals(c,i+1) NURBSderivs(c,i+1)] = ...
            NURBSbasis(i, p, xi(c), uKnot, weights);
    end
end

figure
plot(xi, BsplineVals,'LineWidth',1.4)

%% order elevation to quadratic 

quadraticCurve = nrbdegelev(refineLinearCurve,1); 

% plot basis

uKnot = quadraticCurve.knots;
xi=0:0.001:max(uKnot);
numPts=numel(xi);
m = numel(uKnot)-1;
p = quadraticCurve.order-1;
n = m - p;
weights = ones(1,n);

BsplineVals = zeros(numPts,n);
NURBSderivs = zeros(numPts,n);

for i=1:n
    for c=1:numPts
        [BsplineVals(c,i+1) NURBSderivs(c,i+1)] = ...
            NURBSbasis(i, p, xi(c), uKnot, weights);
    end
end

figure
plot(xi, BsplineVals,'LineWidth',1.4)

%% order elevation to cubic

cubicCurve = nrbdegelev(refineLinearCurve,2); 

% plot basis

uKnot = cubicCurve.knots;
xi=0:0.001:max(uKnot);
numPts=numel(xi);
m = numel(uKnot)-1;
p = cubicCurve.order-1;
n = m - p;
weights = ones(1,n);

BsplineVals = zeros(numPts,n);
NURBSderivs = zeros(numPts,n);

for i=1:n
    for c=1:numPts
        [BsplineVals(c,i+1) NURBSderivs(c,i+1)] = ...
            NURBSbasis(i, p, xi(c), uKnot, weights);
    end
end

figure
plot(xi, BsplineVals,'LineWidth',1.4)

%% order elevation to quartic

quarticCurve = nrbdegelev(refineLinearCurve,3); 

% plot basis

uKnot = quarticCurve.knots;
xi=0:0.001:max(uKnot);
numPts=numel(xi);
m = numel(uKnot)-1;
p = quarticCurve.order-1;
n = m - p;
weights = ones(1,n);

BsplineVals = zeros(numPts,n);
NURBSderivs = zeros(numPts,n);

for i=1:n
    for c=1:numPts
        [BsplineVals(c,i+1) NURBSderivs(c,i+1)] = ...
            NURBSbasis(i, p, xi(c), uKnot, weights);
    end
end

figure
plot(xi, BsplineVals,'LineWidth',1.4)