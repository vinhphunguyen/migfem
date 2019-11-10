% Data file for a center interface crack (bimaterial crack) in an 
% infinite plate problem
% This problem is solved in Sukumar IJNME, 2004.
% Vinh Phu Nguyen
% Hue, Vietnam April 2012

addpath('~/code/xfem-efg-matlab/fem_util');
addpath ../nurbs-geopdes/inst/
addpath ../nurbs-util/
addpath ../meshing/
addpath ../fem-functions/
addpath ../post-processing/
addpath ../xiga/

W = 20;
p = 2;
q = 2;

noPtsX = 11;
noPtsY = 11;

gcoord=meshRectangularCoord(W,W,noPtsX-1,noPtsY-1);
controlPts=gcoord;

weights = ones(1,noPtsX*noPtsY)';

knotUTemp = linspace(0,1,noPtsX-p+1);
knotVTemp = linspace(0,1,noPtsY-q+1);

uKnot = [0 0 knotUTemp 1 1];
vKnot = [0 0 knotVTemp 1 1];

noCtrPts       = noPtsX   * noPtsY;
noDofs         = noCtrPts * 2;

%% generate element connectivity ...

generateIGA2DMesh

% plot the mesh

buildVisualizationMesh;

% Geometry of the material interface

y0 = W/2;
linInterfaceGeo   = [0 y0; W y0];

% Geometry of the crack

noCracks = 1;                 % number of cracks

% data structures for cracks

xCrack   = zeros(noCracks,2,2);
xTips    = zeros(noCracks,2);

a = W/2;                     % crack length
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

incEnrId    = 3; %3;
bimatEnrId  = 4; %4;

enrich_node = zeros(noCtrPts,1);

count1 = 0;

% loop over elements

for iel = 1 : numelem
    sctr    = elementV(iel,:);
    sctrIGA = element(iel,:);
    
    phi     = chi(1,sctr);
        
    if ( max(phi)*min(phi) < 0 )
        count1                 = count1 + 1 ; % ah, one split element
        splitElems(count1)     = iel;                                
        enrich_node(sctrIGA)   = incEnrId;        
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
                enrich_node(sctrIGA) = bimatEnrId; % both tip and interface enr.
                crack_node(sctrIGA)  = iCr;
            end
        end
    end
end

split_nodes = find(enrich_node == 1);
tip_nodes   = find(enrich_node == 2);
itip_nodes  = find(enrich_node == 4);
inc_nodes   = find(enrich_node == 3);
iTip_nodes  = find(enrich_node == 5);

% Plot mesh and enriched nodes to check
figure
hold on
plot_mesh(node,elementV,'Q4','k-');
% plot the interface
plot(linInterfaceGeo(:,1),linInterfaceGeo(:,2),'k-','Linewidth',4.5);
% plot the crack
cr = plot(xCr(:,1),xCr(:,2),'b-');
set(cr,'LineWidth',5);
% plot elements cut by the circle
%plot_mesh(node,elementV(splitElems,:),'Q4','r-');
%plot(gps(:,1),gps(:,2),'+');
n1 = plot(controlPts(inc_nodes,1),controlPts(inc_nodes,2),'r*');
n2 = plot(controlPts(split_nodes,1),controlPts(split_nodes,2),'rs');
n3 = plot(controlPts(itip_nodes,1),controlPts(itip_nodes,2),'ro');
n5 = plot(controlPts(iTip_nodes,1),controlPts(iTip_nodes,2),'rp');
n4 = plot(controlPts(tip_nodes,1),controlPts(tip_nodes,2),'r>');
set(n1,'MarkerSize',16,'LineWidth',1.9);
set(n2,'MarkerSize',16,'LineWidth',1.9);
set(n3,'MarkerSize',16,'LineWidth',1.9);
set(n4,'MarkerSize',16,'LineWidth',1.9);
set(n5,'MarkerSize',16,'LineWidth',1.9);
n6 = plot(controlPts(:,1),controlPts(:,2),'kp');
set(n6,'MarkerSize',16,'LineWidth',1.);