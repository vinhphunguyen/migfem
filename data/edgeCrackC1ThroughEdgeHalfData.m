% data for edge crack problem with C1 elements

p = 1;
q = 1;
L=1;

noPtsX = 30;
noPtsY = 31;

gcoord=meshRectangularCoord(1,1,noPtsX-1,noPtsY-1);
controlPts=gcoord;

weights = ones(1,noPtsX*noPtsY)';

knotUTemp = linspace(0,1,noPtsX-p+1);
knotVTemp = linspace(0,1,noPtsY-q+1);

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

noCracks = 1;                 % number of cracks

% crack data

xCrack   = zeros(noCracks,2,2);
xTips    = zeros(noCracks,2);

D = 0;
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

enrich_node = zeros(noCtrPts,1);
crack_node  = zeros(noCtrPts,1); % which crack to which the node enriched

count1 = 0;
count2 = 0;

tNodes = (p+1)*(p+1);

for i=1:p
    tNodes = [tNodes (p+1)*(p+1)-i]; 
end

tol = 1e-10;

split_elem = [];

% loop over elements

for iel = 1 : numelem
    sctr    = elementV(iel,:);
    sctrIGA = element(iel,:);
    
    % loop over cracks
    for iCr = 1 : noCracks
        phi  = levelSets(iCr,sctr,1);
        psi  = levelSets(iCr,sctr,2);
        
        if ( max(phi)*min(phi) <= tol ) && ( min(phi) < 0)
            if max(psi) < 0
                count1                 = count1 + 1 ; % ah, one split element
                split_elem(count1)     = iel;                
                tip_enr_pos = find(enrich_node(sctrIGA)==2);
                tip_enr     = sctrIGA(tip_enr_pos);
                
                if isempty(tip_enr)
                    enrich_node(sctrIGA(tNodes)) = 1;
                else
                    hea_enr = setdiff(sctrIGA,tip_enr);
                    enrich_node(hea_enr) = 1;
                end
                crack_node(sctrIGA) = iCr;
            end
        end
        
        if ( max(phi)*min(phi) <= tol ) && (max(psi)*min(psi) < 0)
            count2               = count2 + 1 ; % ah, one tip element
            tip_elem(count2)     = iel;
            enrich_node(sctrIGA) = 2;
            crack_node(sctrIGA)  = iCr;
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
            
            