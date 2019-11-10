addpath ../nurbs-geopdes/inst/
addpath ../C_files/
addpath('~/code/xfem-efg-matlab/fem_util');

p       = 3;
knotVec = [0 0 0 0 1 2 3 4 5 5 5 5];
knotVec = 1/max(knotVec)*knotVec;
weights = [1 1 1 1 1 1 1 1 1 1 1 1];

noPts   = 120;
xi      = linspace(0,1,noPts);

BsplineVals = zeros(noPts,length(weights));

for i=1:length(weights)
    for c=1:noPts
        [BsplineVals(c,i+1) NURBSderivs(c,i+1)] = ...
            NURBSbasis(i, p, xi(c), knotVec, weights);
    end
end

figure
hold on
plot(xi,BsplineVals,'Linewidth',2.5)
set(gcf,'color','white')
set(gca,'XTick',[0:1/5:1])

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
exportfig(gcf,'cubicC1Spline.eps',opts)