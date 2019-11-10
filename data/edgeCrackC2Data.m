% data for edge crack problem with C1 elements
% Vinh Phu Nguyen
% Johns Hopkins University

p = 3;
q = 3;
L =1;

noPtsX = 26;
noPtsY = 52;

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

noCracks = 1;                 % number of cracks

xCrack   = zeros(noCracks,2,2);
xTips    = zeros(noCracks,2);

D = 2;
a = 0.45;                     % crack length
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

layerEnrichment = 1;
if (layerEnrichment)  
    x = node(elementV(tip_elem(1),:),:);
    % Area = sum of areas of each sub-triangle
    x0 = x(1,1);
    y0 = x(1,2);
    
    x1 = x(2,1);
    y1 = x(2,2);
    
    x2 = x(3,1);
    y2 = x(3,2);
    
    x3 = x(4,1);
    y3 = x(4,2);
    
    A1 = 0.5 * ((x0-x2)*(y1-y2) - (x1-x2)*(y0-y2)) ;
    A2 = 0.5 * ((x0-x3)*(y2-y3) - (x2-x3)*(y0-y3)) ;
    area = A1 + A2;
    
    % J radius = fac * sqrt(area);
    
    radius = 1.5 * sqrt(area)
    center = xTip;
    
    r=[];
    % Distance from the center of tip element
    for i = 1 : numnode
        sctr = node(i,:);
        rho  = sqrt((sctr(1)-center(1))^2+(sctr(2)-center(2))^2);
        r    = [r,rho];
    end
    test = r-radius;
    test = test(elementV)'; % put nodal test into element test
    % test(4,numelem) for Q4 elements
    test = max(test).*min(test); % test(1,numelem):
    layer = find(test<=0);
    
    for iel = 1 : length(layer)
        e   = layer(iel);
        sctr    = elementV(e,:);
        sctrIGA = element(e,:);
        
        enrich_node(sctrIGA) = 2; % both tip and interface enr.
        crack_node(sctrIGA)  = 1;
    end
    
end

split_nodes = find(enrich_node == 1);
tip_nodes   = find(enrich_node == 2);

% Plot mesh and enriched nodes to check

figure
hold on
plot_mesh(node,elementV,'Q4','b-',1.2);
plot_mesh(node,elementV(layer,:),'Q4','r-',1.2);
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
            
            
