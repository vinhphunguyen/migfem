% data for 3D edge crack problem with C0 elements
% Vinh Phu Nguyen
% Johns Hopkins University

lX = 10;
lY = 2;
lZ = 10;

noPtsX = 16;
noPtsY = 2;
noPtsZ = 16;

% control points come from a structured mesh of trilinear
% brick elements

[controlPts,elementVV]=makeB8mesh(lX,lY,lZ,noPtsX,noPtsY,noPtsZ);

% knot vectors

knotUTemp = linspace(0,1,noPtsX);
knotVTemp = linspace(0,1,noPtsY);
knotWTemp = linspace(0,1,noPtsZ);

uKnot = [0 knotUTemp 1];
vKnot = [0 knotVTemp 1];
wKnot = [0 knotWTemp 1];

% basis order

p = 1;
q = 1;
r = 1;

% weights

weights = ones(1,noPtsX*noPtsY*noPtsZ)';

noCtrPts = noPtsX * noPtsY * noPtsZ;
noDofs   = noCtrPts * 3;

% generate element connectivity ...

generateIGA3DMesh

% crack data

noCracks = 1;                 % number of cracks

% data structures for cracks

a0 = 5;                     % crack length
xCr   = [0 lZ/2; a0 lZ/2];
xTip  = [a0 lZ/2];

% build visualization B8 mesh

buildVisualization3dMesh;

% level set computation

numnode   = size(node,1);
numelem   = size(elementV,1);

% compute level sets on B8 mesh

levelSetsB8 = zeros(numnode,2);

for i=1:numnode
    x = node(i,1); 
    y = node(i,2);
    z = node(i,3);
    
    levelSetsB8(i,1) = z-xCr(1,2);
    levelSetsB8(i,2) = x-a0;
end

% choose enriched nodes

enrich_node = zeros(noCtrPts,1);

count1 = 0;
count2 = 0;

% loop over elements
for iel = 1 : numelem
    sctr    = elementV(iel,:);
    sctrIGA = element(iel,:);
        
    phi  = levelSetsB8(sctr,1);
    psi  = levelSetsB8(sctr,2);
    
    if ( max(phi)*min(phi) < 0 )
        if     max(psi) < 0
            count1 = count1 + 1 ; % ah, one split element
            split_elem(count1)     = iel;
            enrich_node(sctrIGA)   = 1;
        elseif max(psi)*min(psi) < 0
            count2 = count2 + 1 ; % ah, one tip element
            tip_elem(count2)     = iel;
            enrich_node(sctrIGA) = 2;            
        end
    end
    
end

split_nodes = find(enrich_node == 1);
tip_nodes   = find(enrich_node == 2);

% level sets for NURBS mesh

levelSets = levelSetsB8;

%% Plot mesh and enriched nodes to check

figure
hold on
plot_mesh(node,elementV,'B8','b-',1.2);
view(3)
n1 = plot3(controlPts(split_nodes,1),...
           controlPts(split_nodes,2),...
           controlPts(split_nodes,3),'r*');
n2 = plot3(controlPts(tip_nodes,1),...
           controlPts(tip_nodes,2),...
           controlPts(tip_nodes,3),'rs');
       
set(n1,'MarkerSize',16,'LineWidth',1.07);
set(n2,'MarkerSize',16,'LineWidth',1.07);
axis off
set(gcf, 'color', 'white');

% plot(controlPtsX, controlPtsY,'ro',...
%     'MarkerEdgeColor','k',...
%     'MarkerFaceColor','w',...
%     'MarkerSize',9,'LineWidth',1.0);

