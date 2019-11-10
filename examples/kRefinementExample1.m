addpath ../nurbs-geopdes/inst/
addpath ../C_files/
addpath ../fem_util/

% k-refinement= order elevation + knot insertion
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

%% order elevation to quadratic 

quadraticCurve = nrbdegelev(linearCurve,1); 

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

%% knot insertion

newKnots = [1/3 2/3];

refineQuadraticCurve = nrbkntins(quadraticCurve,newKnots);

% plot basis

uKnot = refineQuadraticCurve.knots;
xi=0:0.001:max(uKnot);
numPts=numel(xi);
m = numel(uKnot)-1;
p = refineQuadraticCurve.order-1;
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

cubicCurve = nrbdegelev(linearCurve,2); 

%% knot insertion

newKnots = [1/3 2/3];

refineCubicCurve = nrbkntins(cubicCurve,newKnots);

% plot basis

uKnot = refineCubicCurve.knots;
xi=0:0.001:max(uKnot);
numPts=numel(xi);
m = numel(uKnot)-1;
p = refineCubicCurve.order-1;
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

quarticCurve = nrbdegelev(linearCurve,3); 

%% knot insertion

newKnots = [1/3 2/3];

refineQuarticCurve = nrbkntins(quarticCurve,newKnots);

% plot basis

uKnot = refineQuarticCurve.knots;
xi=0:0.001:max(uKnot);
numPts=numel(xi);
m = numel(uKnot)-1;
p = refineQuarticCurve.order-1;
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

