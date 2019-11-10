addpath ../nurbs-geopdes/inst/
addpath ../C_files/
addpath ../fem_util


uKnot = [0 0 0 0 0.25 0.5 0.75 1 1 1 1];
vKnot = [0 0 0 0 0.25 0.5 0.75 1 1 1 1];
p     = 3;
q     = 3;
weights = [1 1 1 1 1 1 1];

[X,Y] = meshgrid (linspace(0,1,60));
R     = zeros(size(X,1),size(X,1));

for i=1:size(X,1)
    for j=1:size(X,1)
        Xi  = X(1,i);
        Eta = Y(j,1);
        [Ni dN] = NURBSbasis (3, p, Xi,  vKnot, weights);
        [Mi dN] = NURBSbasis (3, q, Eta, uKnot, weights);
        R(i,j) = Ni*Mi;
    end
end

noPts   = 120;
xi      = linspace(0,1,noPts);
N       = zeros(noPts,1);
for i=1:noPts
    [N(i) dN] = NURBSbasis (3, p, xi(i), uKnot, weights);
end

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);

figure
hold on
plot(xi,N,'Linewidth',2.5)

%exportfig(gcf,'cubic-N3',opts)

figure
surf (X,Y,R)
title('Basis function associated to a local knot vector')
hold off

