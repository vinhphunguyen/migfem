%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% debug only
% clear node
% clear elementV
%
% buildVisualizationMesh;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% Build a cracked mesh.
% Intersection of crack and mesh is determined.
% Put duplicate nodes there
% Compute displacement jump at these nodes
% Current status:
%  - only work for linear NURBS (or FEM)
% Vinh Phu Nguyen
% Delft University of Technology/Johns Hopkins University


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add double nodes at intersection points of
% crack with the element (Q4) edges.

noSplitElems = length(split_elem);

newNodeCoord = [];
newNodeIndex = [];

xCr  = reshape(xCrack(1,:,:),2,2);
xTip = xTips(1,:);
y0   = xCr(1,2);

count = 1;

for ie=1:noSplitElems
    elemId = split_elem(ie);
    sctr   = elementV(elemId,:);
    eCoord = node(sctr,:);
    xMin   = min(eCoord(:,1));
    xMax   = max(eCoord(:,1));
    anode  = [xMin y0];
    
    % add two overlapping nodes here
    % they will have different displacements
    
    newNodeCoord = [newNodeCoord; anode];
    newNodeCoord = [newNodeCoord; anode];
    newNodeIndex = [newNodeIndex (numnode + count) (numnode+count+1)];
    count        = count + 2;
end

% tip element

sctr = elementV(tip_elem,:);
sctrTip = sctr; % save the connectivity of tip element
eCoord = node(sctr,:);
xMin   = min(eCoord(:,1));
xMax   = max(eCoord(:,1));
yMin   = min(eCoord(:,2));
yMax   = max(eCoord(:,2));

% add 6 nodes here
%         8
% --------!--------
% |       !      |
% |       !      |
% ========5------|-6
% |       !      |
% |       !      |
% --------!-------
%         7

newNodeCoord = [newNodeCoord; xMin xTip(1,2)];
newNodeCoord = [newNodeCoord; xMin xTip(1,2)];
newNodeCoord = [newNodeCoord; xTip];
newNodeCoord = [newNodeCoord; xMax y0];
newNodeCoord = [newNodeCoord; xTip(1,1) yMin];
newNodeCoord = [newNodeCoord; xTip(1,1) yMax];

id1 = newNodeIndex(end)+1;
id2 = id1+1;
newNodeIndex = [newNodeIndex,id1,id2];

id5 = newNodeIndex(end)+1;
id6 = id5+1;
id7 = id6+1;
id8 = id7+1;

newNodeIndex = [newNodeIndex,id5,id6,id7,id8];

% change the connectivity of tip element
% and add three more elements

elementV(tip_elem,:) = [sctr(1) id7 id5 id5-2];
elementV = [elementV;
    id7 sctr(2) id6 id5;
    id5 id6 sctr(3) id8;
    id5 id8 sctr(4) id5-1];

% add new nodes to the existing nodes

numOldNodes   = length(node);
node          = [node; newNodeCoord];
numTotalNodes = length(node);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute the displacement jumps along the crack at new nodes

jump = [];

for ie=1:noSplitElems
    elemId = split_elem(ie);
    sctr   = elementV(elemId,:);
    sctrIGA= element(elemId,:);
    eCoord = node(sctr,:);
    
    nne    = length(sctrIGA); % no of control pnts/elem
    
    idu    = index(elemId,1);
    idv    = index(elemId,2);
    xiE    = elRangeU(idu,:); % [xi_i,xi_i+1]
    etaE   = elRangeV(idv,:); % [eta_j,eta_j+1]
    
    xMin   = min(eCoord(:,1));
    yMin   = min(eCoord(:,2));
    yMax   = max(eCoord(:,2));
    
    % point in Q4 parent coordinate
    
    eta    = ( 2*y0 - yMin - yMax ) / ( yMax - yMin );
    pt     = [-1 eta];
    
    % corresponding point in NURBS parametric coord
    
    Xi     = parent2ParametricSpace(xiE,pt(1));
    Eta    = parent2ParametricSpace(etaE,pt(2));
    
    [R dRdxi dRdeta] = NURBS2DBasisDers([Xi;Eta],p,q,uKnot,vKnot,weights');
    
    xp    = [xMin y0]-xTip;           % local coordinates
    r     = sqrt(xp(1)*xp(1)+xp(2)*xp(2));
    
    jumpN  = computeJump(R,r,elemId,pos,enrich_node,U);
    jump   = [jump; jumpN];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
% tip element
%%%%%%%%%%%%%%%%%%%%%%%%%%

%sctrTip   = elementV(elemId,:);
eCoord    = node(sctrTip,:);

idu    = index(tip_elem,1);
idv    = index(tip_elem,2);
xiE    = elRangeU(idu,:); % [xi_i,xi_i+1]
etaE   = elRangeV(idv,:); % [eta_j,eta_j+1]

xMin   = min(eCoord(:,1));
yMin   = min(eCoord(:,2));
yMax   = max(eCoord(:,2));

eta    = ( 2*y0 - yMin - yMax ) / ( yMax - yMin );
pt     = [-1 eta];

Xi     = parent2ParametricSpace(xiE,pt(1));
Eta    = parent2ParametricSpace(etaE,pt(2));

[R dRdxi dRdeta] = NURBS2DBasisDers([Xi;Eta],p,q,uKnot,vKnot,weights');

xp    = [xMin y0]-xTip;           % local coordinates
r     = sqrt(xp(1)*xp(1)+xp(2)*xp(2));
jumpN = computeJump(R,r,tip_elem,pos,enrich_node,U);
jump  = [jump; jumpN];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute displacements of new nodes
% along the crack

noNewNode = size(newNodeCoord,1);

for i=1:noSplitElems+1
    h     = 0.5*jump(i,:);
    dispX = [dispX;-h(1)];
    dispY = [dispY;-h(2)];
    dispX = [dispX;h(1)];
    dispY = [dispY;h(2)];
end

% displacements and stresses for new nodes in the TIP element
% do not need to use NURBS interpolation since
% everything is available at nodes of the Q4 mesh.

eCoord    = node(sctrTip,:);
sigmaXXEl = sigmaXX(sctrTip);
sigmaYYEl = sigmaYY(sctrTip);
sigmaXYEl = sigmaXY(sctrTip);
sigmaVMEl = sigmaVM(sctrTip);

xMin   = min(eCoord(:,1));
xMax   = max(eCoord(:,1));
yMin   = min(eCoord(:,2));
yMax   = max(eCoord(:,2));

eta       = ( 2*y0 - yMin - yMax ) / ( yMax - yMin );
pt        = [1 eta];
[N,dNdxi] = lagrange_basis('Q4',pt);
disp6     = N'*[dispX(sctrTip) dispY(sctrTip)];

sigmaXX6  = N'* sigmaXXEl;
sigmaYY6  = N'* sigmaYYEl;
sigmaXY6  = N'* sigmaXYEl;
sigmaVM6  = N'* sigmaVMEl;

xi        = ( 2*xTip(1) - xMin - xMax ) / ( xMax - xMin );
pt        = [xi -1];
[N,dNdxi] = lagrange_basis('Q4',pt);
disp7     = N'*[dispX(sctrTip) dispY(sctrTip)];

sigmaXX7  = N'* sigmaXXEl;
sigmaYY7  = N'* sigmaYYEl;
sigmaXY7  = N'* sigmaXYEl;
sigmaVM7  = N'* sigmaVMEl;

pt        = [xi 1];
[N,dNdxi] = lagrange_basis('Q4',pt);
disp8     = N'*[dispX(sctrTip) dispY(sctrTip)];

sigmaXX8  = N'* sigmaXXEl;
sigmaYY8  = N'* sigmaYYEl;
sigmaXY8  = N'* sigmaXYEl;
sigmaVM8  = N'* sigmaVMEl;

dispX = [dispX;0;disp6(1);disp7(1);disp8(1)];
dispY = [dispY;0;disp6(2);disp7(2);disp8(2)];

% stresses at the crack tip

pt = inverseQ4Mapping (xTip,eCoord); 
[N,dNdxi] = lagrange_basis('Q4',pt);
sigmaXX5  = N'* sigmaXXEl;
sigmaYY5  = N'* sigmaYYEl;
sigmaXY5  = N'* sigmaXYEl;
sigmaVM5  = N'* sigmaVMEl;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% modify connectivity of cut elements
% and add new elements

for ie=1:noSplitElems
    elemId = split_elem(ie);
    sctr   = elementV(elemId,:);
    
    first  = 2*ie-1;
    second = 2*ie;
    third  = first+2;
    fourth = first+3;
    
    elementV(elemId,:) = [sctr(1) sctr(2) newNodeIndex(third) newNodeIndex(first)];
    newElem            = [newNodeIndex(second) newNodeIndex(fourth) sctr(3) sctr(4)];
    
    elementV = [elementV; newElem];
end

%%%
figure
hold on
plot_mesh(node+fac*[dispX dispY],elementV,'Q4','b*-',1);
set(gcf, 'color', 'white');
axis off
%plot_field(node+fac*[dispX dispY],elementV,'Q4',dispY);

oldNodes = [1:numOldNodes];
newNodes = [numOldNodes+1:length(node)];

n1 = plot(node(oldNodes,1)+fac*dispX(oldNodes),...
          node(oldNodes,2)+fac*dispY(oldNodes),'r*');
n2 = plot(node(newNodes,1)+fac*dispX(newNodes),...
          node(newNodes,2)+fac*dispY(newNodes),'rs');
set(n1,'MarkerSize',16,'LineWidth',1.02);
set(n2,'MarkerSize',16,'LineWidth',1.02);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Export to VTU (Paraview)

numAddedNodes = numTotalNodes - numOldNodes;

for i=1:numAddedNodes-4
    sigmaXX = [sigmaXX; 0.];
    sigmaYY = [sigmaYY; 0.];
    sigmaXY = [sigmaXY; 0.];
    sigmaVM = [sigmaVM; 0.];
end

sigmaXX = [sigmaXX; sigmaXX5; sigmaXX6;sigmaXX7;sigmaXX8];
sigmaYY = [sigmaYY; sigmaYY5; sigmaYY6;sigmaYY7;sigmaYY8];
sigmaXY = [sigmaXY; sigmaXY5; sigmaXY6;sigmaXY7;sigmaXY8];
sigmaVM = [sigmaVM; sigmaVM5; sigmaVM6;sigmaVM7;sigmaVM8];

VTKPostProcess(node,elementV,2,'Quad4',vtuCrackFile,...
    [sigmaXX sigmaYY sigmaXY sigmaVM],[dispX dispY]);


