addpath ../nurbs-geopdes/inst/
addpath ../C_files/
addpath('~/code/xfem-efg-matlab/fem_util');
addpath ../fem-functions/

p       = 3;
knotVec = [0 0 0 0 1 2 3 4 5 6 7 8 8 8 8];
weights = ones(length(knotVec)-p-1,1);

noElemsU = length(unique(knotVec))-1;
[elRangeU,elConnU] = buildConnectivity(p,knotVec,noElemsU);

elemId = 4;
connU  = elConnU(elemId,:);

% points to plot basis functions

noPts   = 200;
xi      = linspace(0,8,noPts);
N       = zeros(noPts,length(weights));

for in=1:length(weights)
    for i=1:noPts
        [N(i,in) dN] = NURBSbasis (in, p, xi(i), knotVec, weights);
    end
end

figure
hold on
plot(xi,N,'Linewidth',2.5)
axis([0 8 0 1])
set(gca,'DataAspectRatio',[2 1 1])
set(gca,'YTick',[0:0.5:1])
set(gcf,'color','white')
%axis normal

Nsub = zeros(p+2,noPts);
for i=1:noPts
    c    = 1;
    for j=0:p+1
        [Ni dN]   = NURBSbasis (connU(3), p, 2*xi(i)-j, knotVec, weights);
        Nsub(c,i) = Nsub(c,i) + Ni;
        c         = c + 1;
    end
end

%Nsub(:,1) = 0;

figure
hold on
for i=3:p+2
  plot(xi,Nsub(i,:),'-r','Linewidth',2.2)
end
axis([0 8 0 1])
set(gca,'DataAspectRatio',[2 1 1])
set(gca,'YTick',[0:0.5:1])
set(gcf,'color','white')

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',24);
%exportfig(gcf,'BSplineSubdivision',opts)

Nsub = zeros(p+2,noPts);
for i=1:noPts
    c    = 1;
    for j=0:p+1
        [Ni dN]   = NURBSbasis (9, p, 4*xi(i)-j, knotVec, weights);
        Nsub(c,i) = Nsub(c,i) + Ni;
        c         = c + 1;
    end
end

figure
hold on
for i=1:p+2
  plot(xi,Nsub(i,:),'-r','Linewidth',2.2)
end
axis([0 8 0 1])
set(gca,'DataAspectRatio',[2 1 1])
set(gca,'YTick',[0:0.5:1])
set(gcf,'color','white')