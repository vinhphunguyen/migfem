% data for edge shear crack problem with Ck elements
% h and k-refinement can be used.
% Vinh Phu Nguyen
% Johns Hopkins University

%%

addpath ../nurbs-geopdes/inst/

controlPts          = zeros(4,2,2);

controlPts(1:2,1,1) = [0;0];
controlPts(1:2,2,1) = [7;0];
controlPts(1:2,1,2) = [0;16];
controlPts(1:2,2,2) = [7;16];

controlPts(4,:,:)   = 1;

uKnot = [0 0 1 1];
vKnot = [0 0 1 1];

%% build NURBS object

solid = nrbmak(controlPts,{uKnot vKnot});

% and evaluate order 

solid = nrbdegelev(solid,[3 3]); 

newKnotsX = [0.43 0.55];
newKnotsY = [0.45 0.55];
newKnots  = {newKnotsX newKnotsY};
solid     = nrbkntins(solid,newKnots);

% h-refinement

convert2DNurbs

refineCount = 3;

for i=1:refineCount
    uKnotVectorU = unique(uKnot);
    uKnotVectorV = unique(vKnot);
    
    % new knots along two directions (uniform)
    
    newKnotsX = uKnotVectorU(1:end-1) + 0.5*diff(uKnotVectorU);
    %newKnotsY = uKnotVectorV(1:end-1) + 0.5*diff(uKnotVectorV);
    newKnotsY = [];
    
    for i=1:length(uKnotVectorV)-1
       xi1 = uKnotVectorV(i);
       xi2 = uKnotVectorV(i+1);
       if (xi1-0.5)*(xi2-0.5) < 0 % contains the crack point
           newKnotsY = [newKnotsY xi1-(xi1-xi2)/3 xi1-2*(xi1-xi2)/3];
       else
           newKnotsY = [newKnotsY 0.5*(xi1+xi2)];
       end
    end
    
    newKnots  = {newKnotsX newKnotsY};
    
    % h-refinement
    
    solid     = nrbkntins(solid,newKnots);
    
    uKnot      = cell2mat(solid.knots(1));
    vKnot      = cell2mat(solid.knots(2));
end

convert2DNurbs

noCtrPts       = noPtsX * noPtsY;
noDofs         = noCtrPts * 2;

%% generate element connectivity ...

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
plot_mesh(node,elementV,'Q4','b-',2);
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

