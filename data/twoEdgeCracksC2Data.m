% data for two edge cracks problem with C2 elements
% Vinh Phu Nguyen
% Johns Hopkins University

p = 3;
q = 3;

noPtsX = 28;
noPtsY = 56;

gcoord=meshRectangularCoord(1,2,noPtsX-1,noPtsY-1);
controlPts=gcoord;

weights = ones(1,noPtsX*noPtsY)';

knotUTemp = linspace(0,1,noPtsX-p+1);
knotVTemp = linspace(0,1,noPtsY-q+1);

uKnot = [0 0 0 knotUTemp 1 1 1];
vKnot = [0 0 0 knotVTemp 1 1 1];

noCtrPts       = noPtsX * noPtsY;
noDofs         = noCtrPts * 2;

% sometimes h-refinement process gives NAN new control pts
% simply remove them with the following

controlPts  = controlPts(1:noCtrPts,:);
weights     = weights(1:noCtrPts);

% generate element connectivity ...

generateIGA2DMesh

% crack data

width =1;
D  = 2;
a  = 0.3;                     % crack length
y0 = 1;
y1 = D-y0;

noCracks = 2;                 % number of cracks

xCrack   = zeros(noCracks,2,2);
xTips    = zeros(noCracks,2);

xCrack(1,:,:)   = [0 y0; a y0];
xCrack(2,:,:)   = [1 y1; 1-a y1];

xTips(1,:)    = [a   y0];
xTips(2,:)    = [1-a y1];

% plot the mesh

buildVisualizationMesh;

% level set computation
% Easy implementation: level sets for two cracks at all nodes

numnode   = size(node,1);
numelem   = size(elementV,1);

levelSetCracks

% Choose enriched nodes...

chooseEnrichedNodes

split_nodes = find(enrich_node == 1);
tip_nodes   = find(enrich_node == 2);

% Plot mesh and enriched nodes to check
figure
hold on
plot_mesh(node,elementV,'Q4','b-');
for iCr = 1 : noCracks   
    xCr = reshape(xCrack(iCr,:,:),2,2);
    cr = plot(xCr(:,1),xCr(:,2),'r-');
    set(cr,'LineWidth',3);
end
n1 = plot(controlPts(split_nodes,1),controlPts(split_nodes,2),'r*');
n2 = plot(controlPts(tip_nodes,1),controlPts(tip_nodes,2),'rs');
set(n1,'MarkerSize',16,'LineWidth',1.07);
set(n2,'MarkerSize',16,'LineWidth',1.07);
axis off
set(gcf, 'color', 'white');

plot(controlPtsX, controlPtsY,'ro',...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','w',...
    'MarkerSize',6,'LineWidth',1.0);

