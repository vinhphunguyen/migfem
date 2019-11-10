addpath ../fem_util
addpath ../C_files/
addpath ../data/
addpath ../meshing/
addpath ../post-processing/
addpath ../fem-functions/
addpath ../analytical-solutions/
addpath ../nurbs-geopdes/inst/
addpath ../nurbs-util/

controlPts          = zeros(4,2,2);

controlPts(1:2,1,1) = [0;0];
controlPts(1:2,2,1) = [3;0];
controlPts(1:2,1,2) = [0;3];
controlPts(1:2,2,2) = [3;3];

controlPts(4,:,:)   = 1;

uKnot = [0 0 1 1];
vKnot = [0 0 1 1];

%% build NURBS object

solid = nrbmak(controlPts,{uKnot vKnot});

% and evaluate order 

 solid = nrbdegelev(solid,[1 1]); 

newKnotsX = [0.5];
newKnotsY = [0.5 0.5];
newKnots  = {newKnotsX newKnotsY};
solid     = nrbkntins(solid,newKnots);

%%

convert2DNurbs

noCtrPts       = noPtsX * noPtsY;
noDofs         = noCtrPts * 2;

%% generate element connectivity ...

generateIGA2DMesh
buildVisualizationMesh;

figure
hold on
plot_mesh(node,elementV,'Q4','b-');
plot(controlPts(:,1),controlPts(:,2),'ro',...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',10);
axis off
%%

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
exportfig(gcf,'ret.eps',opts)


