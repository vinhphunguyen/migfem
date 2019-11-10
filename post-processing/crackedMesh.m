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

nodalSupport = zeros(numnode,1);

for in=1:numnode
    nodalSupport(in) = length(find(elementV==in));
end

sigmaXX = zeros(numnode,1);
sigmaYY = zeros(numnode,1);
sigmaXY = zeros(numnode,1);


for e=1:size(elementV,1)
    connect = elementV(e,:);
    for in=1:4
        nid = connect(in);
        sigmaXX(nid) = sigmaXX(nid) + stress(e,in,1);
        sigmaYY(nid) = sigmaYY(nid) + stress(e,in,2);
        sigmaXY(nid) = sigmaXY(nid) + stress(e,in,3);
    end
end

% nodal averaging for stressses

sigmaXX = sigmaXX./nodalSupport;
sigmaYY = sigmaYY./nodalSupport;
sigmaXY = sigmaXY./nodalSupport;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% compute the displacement jumps along the crack at new nodes

jump = [];

for ie=1:noSplitElems
    elemId = split_elem(ie);
    sctr   = elementV(elemId,:);
    eCoord = node(sctr,:);
    
    %elemHeavisideDofs = element_H_dofs(elemId,pos,U);
    elemDisp  = element_disp(elemId,pos,enrich_node,U);
    aDofs  = elemDisp(9:end);
    
    axDofs = aDofs(1:2:8);
    ayDofs = aDofs(2:2:8);
    
    % change position due to different between NURBS
    % element connectivity and Q4 connectivity
    
    e3        = axDofs(3);
    e4        = axDofs(4);
    axDofs(4) = e3;
    axDofs(3) = e4;
    
    e3        = ayDofs(3);
    e4        = ayDofs(4);
    ayDofs(4) = e3;
    ayDofs(3) = e4;
    
    %axDofs = elemHeavisideDofs(1:4);
    %ayDofs = elemHeavisideDofs(5:8);
    
    yMin   = min(eCoord(:,2));
    yMax   = max(eCoord(:,2));
    
    eta     = ( 2*y0 - yMin - yMax ) / ( yMax - yMin );
    pt = [-1 eta];
    [N,dNdxi] = lagrange_basis('Q4',pt);
    
    if ( length(elemDisp) == 2*2*length(sctr) ) % fully H enriched elems
       jump = [jump; 2*N'*[axDofs ayDofs]];
    else % both H and tip enriched element (element ahead tip element)
        axDofs = [aDofs(1) aDofs(3) aDofs(13) aDofs(11)];
        ayDofs = [aDofs(2) aDofs(4) aDofs(14) aDofs(12)];
        jump   = [jump; 2*N'*[axDofs' ayDofs']];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
% tip element
%%%%%%%%%%%%%%%%%%%%%%%%%%

%sctrTip   = elementV(elemId,:);
eCoord    = node(sctrTip,:);
elemDisp  = element_disp(tip_elem,pos,enrich_node,U);
bDofs     = elemDisp(9:end);
bxDofs    = [bDofs(1) bDofs(9) bDofs(17) bDofs(25)];
byDofs    = [bDofs(2) bDofs(10) bDofs(18) bDofs(26)];

xMin   = min(eCoord(:,1));
yMin   = min(eCoord(:,2));
yMax   = max(eCoord(:,2));

eta     = ( 2*y0 - yMin - yMax ) / ( yMax - yMin );
pt = [-1 eta];
[N,dNdxi] = lagrange_basis('Q4',pt);

xp    = [xMin y0]-xTip;           % local coordinates
r     = sqrt(xp(1)*xp(1)+xp(2)*xp(2));

jump = [jump; 2*sqrt(r)*N'*[bxDofs' byDofs']];

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

% displacements and stresses for new nodes
% in the tip element

eCoord    = node(sctrTip,:);
sigmaXXEl = sigmaXX(sctrTip);
sigmaYYEl = sigmaYY(sctrTip);
sigmaXYEl = sigmaXY(sctrTip);

xMin   = min(eCoord(:,1));
xMax   = max(eCoord(:,1));
yMin   = min(eCoord(:,2));
yMax   = max(eCoord(:,2));

eta       = ( 2*y0 - yMin - yMax ) / ( yMax - yMin );
pt        = [1 eta];
[N,dNdxi] = lagrange_basis('Q4',pt);
disp6     = N'*[U(2*sctrTip-1) U(2*sctrTip)];

sigmaXX6  = N'* sigmaXXEl;
sigmaYY6  = N'* sigmaYYEl;
sigmaXY6  = N'* sigmaXYEl;

xi        = ( 2*xTip(1) - xMin - xMax ) / ( xMax - xMin );
pt        = [xi -1];
[N,dNdxi] = lagrange_basis('Q4',pt);
disp7     = N'*[U(2*sctrTip-1) U(2*sctrTip)];

sigmaXX7  = N'* sigmaXXEl;
sigmaYY7  = N'* sigmaYYEl;
sigmaXY7  = N'* sigmaXYEl;

pt        = [xi 1];
[N,dNdxi] = lagrange_basis('Q4',pt);
disp8     = N'*[U(2*sctrTip-1) U(2*sctrTip)];

sigmaXX8  = N'* sigmaXXEl;
sigmaYY8  = N'* sigmaYYEl;
sigmaXY8  = N'* sigmaXYEl;

dispX = [dispX;0;disp6(1);disp7(1);disp8(1)];
dispY = [dispY;0;disp6(2);disp7(2);disp8(2)];

% stresses at the crack tip

pt = inverseQ4Mapping (xTip,eCoord); 
[N,dNdxi] = lagrange_basis('Q4',pt);
sigmaXX5  = N'* sigmaXXEl;
sigmaYY5  = N'* sigmaYYEl;
sigmaXY5  = N'* sigmaXYEl;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Export to VTU (Paraview)

numAddedNodes = numTotalNodes - numOldNodes;

for i=1:numAddedNodes-4
    sigmaXX = [sigmaXX; 0.];
    sigmaYY = [sigmaYY; 0.];
    sigmaXY = [sigmaXY; 0.];
end

sigmaXX = [sigmaXX; sigmaXX5; sigmaXX6;sigmaXX7;sigmaXX8]; 
sigmaYY = [sigmaYY; sigmaYY5; sigmaYY6;sigmaYY7;sigmaYY8]; 
sigmaXY = [sigmaXY; sigmaXY5; sigmaXY6;sigmaXY7;sigmaXY8]; 

VTKPostProcess(node,elementV,2,'Quad4',vtuCrackFile,...
               [sigmaXX sigmaYY sigmaXY],[dispX dispY]);
           
%%%
plot_mesh(node+fac*[dispX dispY],elementV,'Q4','b*-');
plot_field(node+fac*[dispX dispY],elementV,'Q4',dispX);
%plot_field(node+fac*[dispX dispY],elementV,'Q4',sigmaXX);
colorbar
title('Displacement in x direction')
axis off
set(gcf,'color','white')           
