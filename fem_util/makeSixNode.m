% [element,node,edges]=makeSixNode(element,node,edges)
% function makeSixNode: contructs the new element and node matrices 
% for a six node triangular element starting from a three node
% triangular mesh element. In each element three nodes are added.
% element is the connectivity matrix for element
% node is the node matrix
% edges is the list of the connectivity matrix for the boundary
%
% This function must be called after tricheck function
% 
% Northwestern University Feb 24 
% Furio Lorenzo Stazi 


function[element,node,edges]=makeSixNode(element,node,edges)

% debug
% node=[0 0;
%       2 0;
%       0 2;
%       3 3;
%       4 3;
%       4 0];
% element=[1 2 3;
%          3 2 4;
%          4 2 5
%          5 2 6];
% edges={[3 4 
%         4 5];
%        [1 2 
%         2 6]}
% debug

if size(element,1)==1
    onlyOne=1;
else
    onlyOne=0;
end

% create a six node triangular mesh
numnode=size(node,1);
numelem=size(element,1);
newelement=[];
X=node(:,1);
Y=node(:,2);
Xelem=X(element);
Yelem=Y(element);
if size(Xelem,2)==1
    Xelem=Xelem';
    Yelem=Yelem';
end
edge1Nodes=[Xelem(:,2)+Xelem(:,1) Yelem(:,2)+Yelem(:,1)]/2;
edge2Nodes=[Xelem(:,3)+Xelem(:,2) Yelem(:,3)+Yelem(:,2)]/2;
edge3Nodes=[Xelem(:,3)+Xelem(:,1) Yelem(:,3)+Yelem(:,1)]/2;
newConn=(1:numelem)' ;
newConn=[newConn newConn+numelem newConn+2*numelem]+numnode;

node=[node;edge1Nodes;edge2Nodes;edge3Nodes];
element=[element,newConn];

[newnode,inverseMap,newNodeNumberMap]=unique(node,'rows');
node=newnode;
element=newNodeNumberMap(element);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  ADD  NODES IN THE EDGE BOUNDARIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for ed=1:size(edges,1)
    % find new number node for edge nodes
    %edges{ed}=newNodeNumberMap(edges{ed})';
    
    boundary=edges{ed};
    if isempty(boundary)==0
        
        
        Xboun=X(boundary);
        Yboun=Y(boundary);
        if size(boundary,1)==1
            Xboun=Xboun';
            Yboun=Yboun';
        end
        middle=[(Xboun(:,1)+Xboun(:,2))/2 (Yboun(:,1)+Yboun(:,2))/2];
        
        newBoundaryConn=(1:size(boundary,1))'+size(node,1);
        
        node=[node;middle];
        if size(boundary,1)==1
            boundary=[newNodeNumberMap(boundary)' newBoundaryConn];
        else
            boundary=[newNodeNumberMap(boundary) newBoundaryConn];
        end
        
        [newnode,inverseMap,newNodeNumberMap2]=unique(node,'rows');
        node=newnode;
        element=newNodeNumberMap2(element);
        if size(boundary,1)==1
            boundary=newNodeNumberMap2(boundary)';
        else    
            boundary=newNodeNumberMap2(boundary);
        end
        edges{ed}=boundary;
    end
end
    
if onlyOne==1
    element=element';
end

% add the new nodes in the edge boundaries
% loop over the element


% plot debug
% trisurf(element(:,1:3),node(:,1),node(:,2),0*node(:,1))
% hold on
% plot(node(:,1),node(:,2),'ro')
% for ed=1:size(edges,1)
%     hold on
%     plot(node(edges{ed},1),node(edges{ed},2),'b*')
% end
% view(2)
% axis equal
