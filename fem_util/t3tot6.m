function [elements,node]=t3tot6(elements,node)

% [element,node]=makeSixNode(element,node)
%
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

% % test input
% [node,elements]=make_cross_mesh([0 0],[10 10],4,4);
% clf
% plot_mesh(node,elements{5},'T3','g-o')

% create a six node triangular mesh
X=node(:,1);
Y=node(:,2);

% add new nodes to node and element connectivies
for i=1:length(elements)
   
  numnode=size(node,1); 
  element=elements{i};
  numelem=size(element,1);
  
  if ( numelem == 1 )
    Xelem=X(element)';
    Yelem=Y(element)';
  else
    Xelem=X(element);
    Yelem=Y(element);
  end
  
  newConn=(1:numelem)' ;
  if ( size(element,2) == 4 )  % interior element
    
    edge1Nodes=[Xelem(:,2)+Xelem(:,1) Yelem(:,2)+Yelem(:,1)]/2;
    edge2Nodes=[Xelem(:,3)+Xelem(:,2) Yelem(:,3)+Yelem(:,2)]/2;
    edge3Nodes=[Xelem(:,3)+Xelem(:,1) Yelem(:,3)+Yelem(:,1)]/2;
    
    node=[node;edge1Nodes;edge2Nodes;edge3Nodes];
    
    newConn=[newConn newConn+numelem newConn+2*numelem]+numnode;
    
  elseif ( size(element,2) == 3 )  % boundary element
    
    midNodes=[Xelem(:,2)+Xelem(:,1) Yelem(:,2)+Yelem(:,1)]/2;
    
    node=[node;midNodes];
    newConn=[newConn]+numnode;
    
  end
  
  % add nodes to connectivies
  element=[element,newConn];
  elements{i}=element;
  
end

% eliminate duplicate nodes
[newnode,inverseMap,newNodeNumberMap]=unique(node,'rows');
node=newnode;

for i=1:length(elements)
  elements{i}=newNodeNumberMap(elements{i});
  
  if ( size(elements{i},2) == 1 )
    elements{i}=elements{i}';
  end
  
end

