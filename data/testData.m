% Data file for an edge interface crack (bimaterial crack) problem
% This problem is solved in Karihaloo IJNME, 2004.

addpath ../fem_util/;
addpath ../nurbs-geopdes/inst/
addpath ../nurbs-util/
addpath ../meshing/
addpath ../fem-functions/
addpath ../post-processing/
addpath ../xiga/
addpath ../nurbs-geopdes/inst/

L = 10;

p = 1;
q = 1;

pp = 2; % raised order
qq = 1;

noPtsX = 5;
noPtsY = 5;

gcoord=meshRectangularCoord(L,L,noPtsX-1,noPtsY-1);
controlPts=gcoord;

weights = ones(1,noPtsX*noPtsY)';

knotUTemp = linspace(0,1,noPtsX);
knotVTemp = linspace(0,1,noPtsY);

uKnot = [0 knotUTemp 1];
vKnot = [0 knotVTemp 1];

%% convert to NURBS object for knot insertion

controlPts1          = zeros(4,noPtsX,noPtsY);
for i=1:noPtsY
    en  = noPtsX*(i-1)+1;
    st  = noPtsX*i;
    controlPts1(1:2,1:noPtsX,i) = controlPts(en:st,:)';
end
controlPts1(4,:,:)   = 1;

solid = nrbmak(controlPts1,{uKnot vKnot});

% and evaluate order 

solid = nrbdegelev(solid,[pp qq]); 

%% convert to IGA format

convert2DNurbs
noCtrPts       = noPtsX * noPtsY;
noDofs         = noCtrPts * 2;
    

%% generate element connectivity ...

generateIGA2DMesh

% plot the mesh

buildVisualizationMesh;

% Plot mesh and enriched nodes to check
figure
hold on
plot_mesh(node,elementV,'Q4','b-');
n5 = plot(controlPts(:,1),controlPts(:,2),'r*');
set(n5,'MarkerSize',12,'LineWidth',1.01);

%% write the mesh to file

fileName  = 'test.nurbs';
file      = fopen(fileName, 'wt');

%% Write headers
fprintf(file, 'NumberOfNodess %g NumberOfElements %g \n', noCtrPts, noElems);

%% Write point data
fprintf(file, 'NODES \n');

for i=1:noCtrPts
    fprintf(file, ' %g %f ',  i, controlPts(i,1:2));
    fprintf(file, '\n');
end

%% Print cell connectivity
fprintf(file, 'ELEMENTS \n');

for i=1:noElems
    fprintf(file, '%g ',  element(i,:));
    fprintf(file, '\n');
end

%% Write material id

fprintf(file, 'MATERIAL ID\n');

matID = 1;

for i=1:noElems
    if i > (noElemsU * noElemsV/2) 
        matID = 2; 
    end
    fprintf(file, '%g %g',  i, matID);
    fprintf(file, '\n');
end

bottomNodes = find(controlPts(:,2)==0)';
rightNodes  = find(controlPts(:,1)==L)';
leftNodes   = find(controlPts(:,1)==0)';
topNodes    = find(controlPts(:,2)==L)';

%% Write node groups
fprintf(file, 'NODE GROUPS %g \n', 3);
fprintf(file, 'bottomNodes\n');
for i=1:length(bottomNodes)
  fprintf(file, '%g ', bottomNodes(i));
end
fprintf(file,'\n');

fprintf(file, 'topNodes\n');
for i=1:length(topNodes)
  fprintf(file, '%g ', topNodes(i));
end
fprintf(file,'\n');

fprintf(file, 'rightNodes\n');
for i=1:length(rightNodes)
  fprintf(file, '%g ', rightNodes(i));
end
fprintf(file,'\n');
