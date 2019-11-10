addpath ../nurbs-geopdes/inst/
addpath ../C_files/
addpath ../fem_util

% plot Bernstein polynomials in [-1,1] for orders from 1 to 3.
% VP Nguyen

clear all

p       = 1;

noPts   = 50;
xi      = linspace(-1,1,noPts);

for i=1:noPts
    for ii=1:p+1
        coef    = nchoosek(p,ii-1);
        N(ii,i) = 1/2^p*coef*(1-xi(i))^(p+1-ii)*(1+xi(i))^(ii-1);
    end
end

figure
hold on
plot(xi,N,'Linewidth',2.5)
set(gca,'XTick',[-1:0.4:1])
set(gca,'YTick',[0:0.2:1])

%%%%%%%%%

clear N

p       = 2;

for i=1:noPts
    for ii=1:p+1
        coef    = nchoosek(p,ii-1);
        N(ii,i) = 1/2^p*coef*(1-xi(i))^(p+1-ii)*(1+xi(i))^(ii-1);
    end
end

figure
hold on
plot(xi,N,'Linewidth',2.5)
set(gca,'XTick',[-1:0.4:1])
set(gca,'YTick',[0:0.2:1])


%%%%%%%%%

clear N

p       = 3;

for i=1:noPts
    for ii=1:p+1
        coef    = nchoosek(p,ii-1);
        N(ii,i) = 1/2^p*coef*(1-xi(i))^(p+1-ii)*(1+xi(i))^(ii-1);
    end
end

figure
hold on
plot(xi,N,'Linewidth',2.5)
set(gca,'XTick',[-1:0.4:1])
set(gca,'YTick',[0:0.2:1])


%%%%%%%%%

clear N

p       = 4;

for i=1:noPts
    for ii=1:p+1
        coef    = nchoosek(p,ii-1);
        N(ii,i) = 1/2^p*coef*(1-xi(i))^(p+1-ii)*(1+xi(i))^(ii-1);
    end
end

figure
hold on
plot(xi,N,'Linewidth',2.5)
set(gca,'XTick',[-1:0.4:1])
set(gca,'YTick',[0:0.2:1])

%plot(xi,Nsub(1,:),'-r','Linewidth',1.0)


%legend('Original','Subdivision')

[X,Y] = meshgrid (linspace(-1,1,60));
R     = zeros(size(X,1),size(X,1));
ii = 2;
jj = 1;
p  = 3;
q  = 1;

for i=1:size(X,1)
    for j=1:size(X,1)
        Xi  = X(1,i);
        Eta = Y(j,1);
        Ni  = 1/2^p*coef*(1-Xi)^(p+1-ii)*(1+Xi)^(ii-1);
        Mi  = 1/2^q*coef*(1-Eta)^(q+1-jj)*(1+Eta)^(jj-1);
        R(i,j) = Ni*Mi;
    end
end

figure
surf (X,Y,R)

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',24);
exportfig(gcf,'BSplineSubdivision.eps',opts)


