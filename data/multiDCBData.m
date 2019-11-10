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

w = 0.1325;
L = 120;

%%

noLayer = 1;

p = 1;
q = 1;

pp = 2;
qq = 1;

noPtsX = 241;
noPtsY = 2;

noPtsX0 = noPtsX;

gcoord1=meshRectangularCoord(L,10*w,noPtsX-1,noPtsY-1);
gcoord2=meshRectangularCoord(L,2*w, noPtsX-1,1);
gcoord3=meshRectangularCoord(L,12*w,noPtsX-1,noPtsY-1);

gcoord2(:,2) = gcoord2(:,2) + 10*w;
gcoord3(:,2) = gcoord3(:,2) + 12*w;


controlPts=[gcoord1;gcoord2(noPtsX+1:end,:);gcoord3(noPtsX+1:end,:)];

noPtsY = 4;

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

fileName  = 'multiDCB.nurbs';
file      = fopen(fileName, 'wt');
file1     = fopen('dcb.mesh', 'wt');

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
    if i <= (noLayer*(noPtsX0-1)) 
        matID = 1; 
        fprintf(file, '%g %g',  i, matID);
        fprintf(file, '\n');
        continue
    end
    
    if i > (noLayer*(noPtsX0-1)) && i <= ((noLayer+1)*(noPtsX0-1))
        matID = 2; 
        fprintf(file, '%g %g',  i, matID);
        fprintf(file, '\n');
        continue
    end
    
    if i > (noLayer+1)*(noPtsX0-1)
        matID = 3;   
        fprintf(file, '%g %g',  i, matID);
        fprintf(file, '\n');
        continue
    end    
end

leftNodes   = find(controlPts(:,1)==0)';
rightNodes  = find(controlPts(:,1)==L)';

%% Write node groups
fprintf(file, 'NODE GROUPS %g \n', 2);
fprintf(file, 'rightNodes\n');
for i=1:length(rightNodes)
  fprintf(file, '%g ', rightNodes(i));
end
fprintf(file,'\n');

fprintf(file, 'force1\n');
fprintf(file, '%g', leftNodes(1));

fprintf(file,'\n');

fprintf(file, 'force2\n');
fprintf(file, '%g', leftNodes(end));

fprintf(file,'\n');
    


