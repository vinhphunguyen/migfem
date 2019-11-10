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

clear all

L = 100;
w = 1.5;

p = 1;
q = 1;

pp = 3;
qq = 0;

noPtsX = 71;
noPtsY = 3;

gcoord=meshRectangularCoord(L,2*w,noPtsX-1,noPtsY-1);
controlPts=gcoord;
controlPts(:,2) = controlPts(:,2) - w;

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

% controlPts          = zeros(4,2,2);
% 
% controlPts(1:2,1,1) = [0;-w];
% controlPts(1:2,2,1) = [L;-w];
% controlPts(1:2,1,2) = [0;w];
% controlPts(1:2,2,2) = [L;w];
% 
% controlPts(4,:,:)   = 1;
% 
% uKnot = [0 0 1 1];
% vKnot = [0 0 1 1];
% 
% %% build NURBS object
% 
% solid = nrbmak(controlPts,{uKnot vKnot});
% 
% % h-refinement
% 
% refineCountX = 6;
% refineCountY = 1;
% newKnotsX = [];
% 
% for i=1:refineCountY
%     uKnotVectorU = unique(uKnot);
%     uKnotVectorV = unique(vKnot);
%    
%     newKnotsY = [];
%     
%     % new knots along two directions (uniform)
%     
%     %newKnotsX = uKnotVectorU(1:end-1) + 0.5*diff(uKnotVectorU);
%     newKnotsY = uKnotVectorV(1:end-1) + 0.5*diff(uKnotVectorV);
%     
%     newKnots  = {newKnotsX newKnotsY};
%     
%     % h-refinement
%     
%     solid     = nrbkntins(solid,newKnots);
%     
%     uKnot      = cell2mat(solid.knots(1));
%     vKnot      = cell2mat(solid.knots(2));
% end
% 
% for i=1:refineCountX
%     uKnotVectorU = unique(uKnot);
%     uKnotVectorV = unique(vKnot);
%    
%     newKnotsY = [];
%     
%     % new knots along two directions (uniform)
%     
%     %newKnotsX = uKnotVectorU(1:end-1) + 0.5*diff(uKnotVectorU);
%     newKnotsX = uKnotVectorU(1:end-1) + 0.5*diff(uKnotVectorU);
%     
%     newKnots  = {newKnotsX newKnotsY};
%     
%     % h-refinement
%     
%     solid     = nrbkntins(solid,newKnots);
%     
%     uKnot      = cell2mat(solid.knots(1));
%     vKnot      = cell2mat(solid.knots(2));
% end
% 
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

fileName  = 'mmb.nurbs';
file      = fopen(fileName, 'wt');
file1     = fopen('mmb.mesh', 'wt');

%% Write headers
fprintf(file, 'NumberOfNodess %g NumberOfElements %g \n', noCtrPts, noElems);

%% Write point data
fprintf(file, 'NODES \n');
fprintf(file1, '<Nodes> \n');

for i=1:noCtrPts
    fprintf(file, ' %g %f ',  i, controlPts(i,1:2));
    fprintf(file, '\n');
    
    fprintf(file1, ' %g %f ',  i, controlPts(i,1:2));
    fprintf(file1, ';\n');
end

fprintf(file1, '</Nodes> \n');

%% Print cell connectivity
fprintf(file, 'ELEMENTS \n');

fprintf(file1, '<Elements> \n');

for i=1:noElems
    fprintf(file, '%g ',  element(i,:));
    fprintf(file, '\n');
    
    fprintf(file1, '%g %g ', i,  element(i,:));
    fprintf(file1, ';\n');
end

fprintf(file1, '</Elements> \n');

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

bottomNodes = find(controlPts(:,2)==-w)';
topNodes    = find(controlPts(:,2)==w)';


%% Write node groups

id = round(length(topNodes)/2);

fprintf(file, 'NODE GROUPS %g \n', 2);
fprintf(file, 'firstNodes\n');
fprintf(file, '%g ', bottomNodes(1));
fprintf(file,'\n');
fprintf(file, 'secondNodes\n');
fprintf(file, '%g ', bottomNodes(end));
fprintf(file,'\n');
fprintf(file, 'force2Nodes\n');
fprintf(file, '%g ', topNodes(id));
fprintf(file,'\n');
fprintf(file, 'force1Nodes\n');
fprintf(file, '%g ', topNodes(end));
fprintf(file,'\n');

    


