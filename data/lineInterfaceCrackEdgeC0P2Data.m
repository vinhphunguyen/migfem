% Data file for an edge interface crack (bimaterial crack) problem
% This problem is solved in Karihaloo IJNME, 2004.

addpath('~/code/xfem-efg-matlab/fem_util');
addpath ../nurbs-geopdes/inst/
addpath ../nurbs-util/
addpath ../meshing/
addpath ../fem-functions/
addpath ../post-processing/
addpath ../xiga/
addpath ../nurbs-geopdes/inst/

L = 3;
p = 2;
q = 2;

noPtsX = 31;
noPtsY = 93;

gcoord     = meshRectangularCoord(L,3*L,noPtsX-1,noPtsY-1);
controlPts = gcoord;

weights    = ones(1,noPtsX*noPtsY)';

knotUTemp = linspace(0,1,noPtsX-p+1);
knotVTemp = linspace(0,1,noPtsY-q+1);

uKnot = [0 0 knotUTemp 1 1];
vKnot = [0 0 knotVTemp 1 1];


%% convert to NURBS object for knot insertion

controlPts1          = zeros(4,noPtsX,noPtsY);
for i=1:noPtsY
    en  = noPtsX*(i-1)+1;
    st  = noPtsX*i;
    controlPts1(1:2,1:noPtsX,i) = controlPts(en:st,:)';
end
controlPts1(4,:,:)   = 1;

solid = nrbmak(controlPts1,{uKnot vKnot});


newKnotsX = knotUTemp(2:end-1);  
newKnotsY = knotVTemp(2:end-1);  

newKnots  = {newKnotsX newKnotsY};

% h-refinement

solid     = nrbkntins(solid,newKnots);

%% convert to IGA format

convert2DNurbs
p = 2;
q = 2;
noCtrPts       = noPtsX * noPtsY;
noDofs         = noCtrPts * 2;
    

%% generate element connectivity ...

generateIGA2DMesh

% plot the mesh

buildVisualizationMesh;

% Geometry of the material interface

y0 = 1.5*L;
linInterfaceGeo   = [0 y0; L y0];

% Geometry of the crack

noCracks = 1;                 % number of cracks

% data structures for cracks

xCrack   = zeros(noCracks,2,2);
xTips    = zeros(noCracks,2);

a = 0.7*L;                     % crack length
xCr   = [0 y0; a y0];
xTip  = [a y0];
seg   = xCr(2,:) - xCr(1,:);   % tip segment
alpha = atan2(seg(2),seg(1));  % inclination angle
QT    =[cos(alpha) sin(alpha); -sin(alpha) cos(alpha)];

xCrack(1,:,:) = xCr;
xTips(1,:)    = xTip; 

% Compute level sets, for interface

numnode = size(node,1);
numelem = size(elementV,1);
nCPts   = size(controlPts,1);
chi     = zeros(1,numnode);
CHI     = zeros(1,nCPts);

x0   = linInterfaceGeo(1,1); y0 = linInterfaceGeo(1,2);
x1   = linInterfaceGeo(2,1); y1 = linInterfaceGeo(2,2);

for i = 1 : numnode
    x                  = node(i,1);
    y                  = node(i,2);
    l                  = sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)) ;
    phi                = (y0-y1)*x + (x1-x0)*y + (x0*y1-x1*y0);
    chi(1,i) = phi/l;
end

for i = 1 : nCPts
    x                  = controlPts(i,1);
    y                  = controlPts(i,2);
    l                  = sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)) ;
    phi                = (y0-y1)*x + (x1-x0)*y + (x0*y1-x1*y0);
    CHI(1,i) = phi/l;
end

% compute level sets for cracks

levelSetCracks

% choose enriched nodes (for both material interface
% and bi mat crack

enrich_node = zeros(noCtrPts,1);

count1 = 0;

cornerNodes = [1 3 7 9];
%cornerNodes = [1:size(element,2)];

% loop over elements

for iel = 1 : numelem
    sctr    = elementV(iel,:);
    sctrIGA = element(iel,:);
    
    phi  = chi(1,sctr);
        
    if ( max(phi)*min(phi) < 0 )
        count1                 = count1 + 1 ; % ah, one split element
        splitElems(count1)     = iel;                                
        enrich_node(sctrIGA(cornerNodes))   = inclusionEnrId;        
    end
end

count1 = 0;
count2 = 0;

% loop over elements
for iel = 1 : numelem
    sctr    = elementV(iel,:);
    sctrIGA = element(iel,:);
    
    % loop over cracks
    for iCr = 1 : noCracks
        phi  = levelSets(iCr,sctr,1);
        psi  = levelSets(iCr,sctr,2);
        
        if ( max(phi)*min(phi) < 0 )
            if max(psi) < 0
                count1                 = count1 + 1 ; % ah, one split element
                split_elem(count1)     = iel;                
                tip_enr_pos = find(enrich_node(sctrIGA)==2);
                tip_enr     = sctrIGA(tip_enr_pos);
                
                if isempty(tip_enr)
                    enrich_node(sctrIGA) = 1;
                else
                    hea_enr = setdiff(sctrIGA,tip_enr);
                    enrich_node(hea_enr) = 1;
                end
                crack_node(sctrIGA) = iCr;
            elseif max(psi)*min(psi) < 0
                count2               = count2 + 1 ; % ah, one tip element
                tip_elem(count2)     = iel;
                enrich_node(sctrIGA) = tipEnrId; % both tip and interface enr.
                crack_node (sctrIGA) = iCr;
            end
        end
    end
end

if (layerEnrichment)  
    layer=[tip_elem-1; tip_elem+1;tip_elem-noElemsU;tip_elem+(noElemsU);...
        tip_elem-1-(noElemsU);tip_elem-1+(noElemsU);...
        tip_elem+1-(noElemsU);tip_elem+1+(noElemsU)];
    
    for iel = 1 : length(layer)
        e   = layer(iel);
        sctr    = elementV(e,:);
        sctrIGA = element(e,:);
        
        enrich_node(sctrIGA) = tipEnrId; % both tip and interface enr.
        crack_node(sctrIGA)  = 1;
    end
end

split_nodes = find(enrich_node == 1);
tip_nodes   = find(enrich_node == 2);
itip_nodes  = find(enrich_node == 4);
iTip_nodes  = find(enrich_node == 5);
inc_nodes   = find(enrich_node == 3);
iTIp_nodes  = find(enrich_node == 6);

% Plot mesh and enriched nodes to check
figure
hold on
plot_mesh(node,elementV,'Q4','b-');
% plot the interface
plot(linInterfaceGeo(:,1),linInterfaceGeo(:,2),'k-','Linewidth',1.9);
% plot the crack
cr = plot(xCr(:,1),xCr(:,2),'cy-');
set(cr,'LineWidth',5);
% plot elements cut by the circle
%plot_mesh(node,elementV(splitElems,:),'Q4','r-');
%plot(gps(:,1),gps(:,2),'+');
n1 = plot(controlPts(inc_nodes,1),controlPts(inc_nodes,2),'r>');
n2 = plot(controlPts(split_nodes,1),controlPts(split_nodes,2),'rs');
n3 = plot(controlPts(itip_nodes,1),controlPts(itip_nodes,2),'ro');
n4 = plot(controlPts(iTip_nodes,1),controlPts(iTip_nodes,2),'rp');
n6 = plot(controlPts(iTIp_nodes,1),controlPts(iTIp_nodes,2),'rp');
n5 = plot(controlPts(:,1),controlPts(:,2),'r*');
set(n1,'MarkerSize',16,'LineWidth',1.01);
set(n2,'MarkerSize',16,'LineWidth',1.01);
set(n3,'MarkerSize',16,'LineWidth',1.01);
set(n4,'MarkerSize',16,'LineWidth',1.01);
set(n6,'MarkerSize',16,'LineWidth',1.01);
set(n5,'MarkerSize',12,'LineWidth',1.01);

c0line = 1;

