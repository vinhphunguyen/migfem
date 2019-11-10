% data for center crack problem with C0 elements
% This problem involves two crack tips.
% NURBS data: control points, knots, basis orders and weights
% Crack data, generate the NURBS mesh, the visualization Q4 mesh
% then level sets computation, and finally enrichment detection
% Vinh Phu Nguyen, Johns Hopkins University

p = 1;
q = 1;

noPtsX = 40;
noPtsY = 40*3;

gcoord=meshRectangularCoord(2,6,noPtsX-1,noPtsY-1);
controlPts=gcoord;

weights = ones(1,noPtsX*noPtsY)';

knotUTemp = linspace(0,1,noPtsX);
knotVTemp = linspace(0,1,noPtsY);

uKnot = [0 knotUTemp 1];
vKnot = [0 knotVTemp 1];

noCtrPts       = noPtsX * noPtsY;
noDofs         = noCtrPts * 2;

% sometimes h-refinement process gives NAN new control pts
% simply remove them with the following 

controlPts  = controlPts(1:noCtrPts,:);
weights     = weights(1:noCtrPts);

% generate element connectivity ...

generateIGA2DMesh

% crack data

D = 6;
a = 0.25;                     % half crack length
xCr   = [0.75 D/2; 1.25 D/2];
xTip  = [0.75 D/2;
         1.25 D/2];
seg1   = xCr(1,:) - xCr(2,:);   % tip segment 1
seg2   = xCr(2,:) - xCr(1,:);   % tip segment 2

alpha1 = atan2(seg1(2),seg1(1));  % inclination angle
alpha2 = atan2(seg2(2),seg2(1));  % inclination angle

%QT    =[cos(alpha) sin(alpha); -sin(alpha) cos(alpha)];


% plot the mesh

buildVisualizationMesh;

% level set computation

x0  = xCr(1,1); y0 = xCr(1,2);
x1  = xCr(2,1); y1 = xCr(2,2);
t1  = 1/norm(seg1)*seg1;
t2  = 1/norm(seg2)*seg2;

numnode = size(node,1);
numelem = size(elementV,1);

for i = 1 : numnode
    x = node(i,1);
    y = node(i,2);
    l   = sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)) ;
    phi = (y0-y1)*x + (x1-x0)*y + (x0*y1-x1*y0);
    ls(i,1) = phi/l;            % normal LS
    ls(i,2) = ([x y]-xTip(1,:))*t1';  % tangent LS tip 1
    ls(i,3) = ([x y]-xTip(2,:))*t2';  % tangent LS tip 1
end

% Choose enriched nodes...

% for one element, if max(phi)*min(phi) < 0
% and max(psi) < 0, then it is a split element
% If max(phi)*min(phi) < 0 and max(psi)*min(psi) < 0, it is
% tip element

% Data structures for elements cut by crack
% Array split_elem contains the number of elements which are completely
% cut by crack. Similarly, array tip_elem stores number of tip element

enrich_node = zeros(noCtrPts,1);
tip_node    = zeros(noCtrPts,1);

count1 = 0;
count2 = 0;

% first Heaviside enriched nodes

for iel = 1 : numelem
    sctr    = elementV(iel,:);
    sctrIGA = element(iel,:);
    phi  = ls(sctr,1);
    psi1 = ls(sctr,2);
    psi2 = ls(sctr,3);
    if ( max(phi)*min(phi) < 0 ) % all elements cut by extended crack
        if max(psi1) < 0 & max(psi2) < 0
            count1 = count1 + 1 ; % ah, one split element
            split_elem(count1) = iel;
            enrich_node(sctrIGA)   = 1;           
        end
    end
end

% then tip enriched nodes, otherwise, some tip enriched
% nodes can be overwritten by H enriched ones

for iel = 1 : numelem
    sctr    = elementV(iel,:);
    sctrIGA = element(iel,:);
    phi  = ls(sctr,1);
    psi1 = ls(sctr,2);
    psi2 = ls(sctr,3);
    if ( max(phi)*min(phi) < 0 ) % all elements cut by extended crack
        if max(psi1)*min(psi1) < 0
            count2 = count2 + 1 ; % ah, one tip 1 element
            tip_elem(count2)      = iel;
            enrich_node(sctrIGA)  = 2;
            tip_node(sctrIGA)     = 1;
        elseif max(psi2)*min(psi2) < 0
            count2 = count2 + 1 ; % ah, one tip 2 element
            tip_elem(count2)     = iel;
            enrich_node(sctrIGA) = 2;            
            tip_node(sctrIGA)    = 2;
        end
    end
end

split_nodes = find(enrich_node == 1);
tip_nodes   = find(enrich_node == 2);

% Plot mesh and enriched nodes to check
figure
hold on
plot_mesh(node,elementV,'Q4','b-');
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
                'MarkerSize',9,'LineWidth',1.0);
            
            
