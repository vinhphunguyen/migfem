addpath ../nurbs-geopdes/inst/
addpath ../C_files/
addpath ../fem_util/

p       = 3;
knotVec = [0 1 2 3 4 5 6 7 8];
weights = [1 1 1 1 1];

noPts   = 120;
xi      = linspace(0,5,noPts);

for i=1:noPts
    [N(i) dN] = NURBSbasis (1, p, xi(i), knotVec, weights);
end
N(1)=0;

Nsub = zeros(p+2,noPts);
for i=1:noPts
    c    = 1;
    for j=0:p+1
        [Ni dN]   = NURBSbasis (1, p, 2*xi(i)-j, knotVec, weights);
        %[Ni dN]   = NURBSbasis (2*1+j, p, 2*xi(i)-j, knotVec, weights);
        coef      = nchoosek(p+1,j);
        Nsub(c,i) = Nsub(c,i) + coef*Ni;
        c         = c + 1;
    end
end

Nsub = 2^(-p)*Nsub;
Nsub(:,1) = 0;

figure
hold on
plot(xi,N,'Linewidth',2.5)

plot(xi,Nsub(1,:),'-r','Linewidth',1.0)
plot(xi,Nsub(2,:),'-r','Linewidth',2.0)
plot(xi,Nsub(3,:),'-r','Linewidth',1.0)
plot(xi,Nsub(4,:),'-r','Linewidth',2.0)
plot(xi,Nsub(5,:),'-r','Linewidth',1.0)

legend('Original','Subdivision')
set(gca,'XTick',[0:1:4])

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',24);
exportfig(gcf,'BSplineSubdivision.eps',opts)

