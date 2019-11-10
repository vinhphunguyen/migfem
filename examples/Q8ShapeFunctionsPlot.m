addpath ../nurbs-geopdes/inst/
addpath ../C_files/
addpath ../fem_util



[X,Y] = meshgrid (linspace(-1,1,100));
N1     = zeros(size(X,1),size(X,1));
N5     = zeros(size(X,1),size(X,1));

for i=1:size(X,1)
    for j=1:size(X,1)
        Xi  = X(1,i);
        Eta = Y(j,1);
        N1(i,j) = 0.25*(1-Xi)*(1-Eta)*(-1-Xi-Eta);
        N5(i,j) = 0.5*(1-Xi*Xi)*(1-Eta);
    end
end


opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);


%%
figure
h=surf (X,Y,N1)
%title('Basis function associated to a local knot vector')
hold off
view(0,90)
colormap(hsv(8))
colorbar
axis equal
grid off
set(gca,'visible','off')
%set(h,'facecolor',[0 .6 0],'linestyle','none')
exportfig(gcf,'Q8-N1.eps',opts)
%%
figure
surf (X,Y,N5)
%title('Basis function associated to a local knot vector')
hold off
view(0,90)
colormap(hsv(8))
colorbar
axis equal
grid off
set(gca,'visible','off')
exportfig(gcf,'Q8-N5.eps',opts)
