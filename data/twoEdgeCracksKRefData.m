% data for two edge cracks problem with Ck elements
% h and k-refinement can be used.
% Vinh Phu Nguyen
% Johns Hopkins University

%%

addpath nurbs-geopdes/inst/

controlPts          = zeros(4,2,2);

controlPts(1:2,1,1) = [0;0];
controlPts(1:2,2,1) = [1;0];
controlPts(1:2,1,2) = [0;2];
controlPts(1:2,2,2) = [1;2];

controlPts(4,:,:)   = 1;

uKnot = [0 0 1 1];
vKnot = [0 0 1 1];

%% build NURBS object

solid = nrbmak(controlPts,{uKnot vKnot});

% and evaluate order 

solid = nrbdegelev(solid,[2 2]); 

newKnotsX = [0.25 0.5 0.75];
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
plot_mesh(node,elementV,'Q4','b-',1);
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

