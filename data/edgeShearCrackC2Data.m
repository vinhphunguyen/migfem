% data for edge crack problem with C2 elements

p = 3;
q = 3;

noPtsX = 28;
noPtsY = 48;

gcoord=meshRectangularCoord(7,16,noPtsX-1,noPtsY-1);
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

noCracks = 1;                 % number of cracks

% crack data

xCrack   = zeros(noCracks,2,2);
xTips    = zeros(noCracks,2);

D = 16;
a = 3.5;                     % crack length
xCr   = [0 D/2; a D/2];
xTip  = [a D/2];
seg   = xCr(2,:) - xCr(1,:);   % tip segment
alpha = atan2(seg(2),seg(1));  % inclination angle
QT    =[cos(alpha) sin(alpha); -sin(alpha) cos(alpha)];

xCrack(1,:,:) = xCr;
xTips(1,:)    = xTip; 


% plot the mesh

buildVisualizationMesh;

% level set computation

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
plot_mesh(node,elementV,'Q4','b-',2);
cr = plot(xCr(:,1),xCr(:,2),'r-');
set(cr,'LineWidth',3);
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
            
            