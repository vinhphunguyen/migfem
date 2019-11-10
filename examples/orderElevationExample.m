addpath ../nurbs-geopdes/inst/
addpath ../C_files/

%% original curve (quadratic)

controlPts = zeros(2,3);

controlPts = [0 0.5 1;
              0 1 0];
          
uKnot      = [0 0 0 1 1 1];          

curve      = nrbmak(controlPts, uKnot);

hold on

plot(controlPts(1,:),controlPts(2,:),'r-o',...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',10);
            
nrbplot(curve,40); 
axis([0 1 0 1])

% plot basis

xi=0:0.01:max(uKnot);
numPts=numel(xi);
m = numel(uKnot)-1;
p = curve.order-1;
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

%% elevated curve

ecurve = nrbdegelev(curve,1); 

figure
hold on

eCtrPoints = ecurve.coefs(1:2,:);

plot(eCtrPoints(1,:),eCtrPoints(2,:),'r-o',...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',10);
            
nrbplot(ecurve,40); 
axis([0 1 0 1])
% plot basis

uKnot = ecurve.knots;
xi=0:0.01:max(uKnot);
numPts=numel(xi);
m = numel(uKnot)-1;
p = ecurve.order-1;
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