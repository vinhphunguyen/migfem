function [node,elemLst]=remove_free_nodes(node,elemLst)
  
%  [node,elemLst]=remove_free_nodes(node,elemLst)
%  
%  This function looks for free nodes and removes then renumbers the nodes
%  in NODE so that the numbering is continous.  Also teh connectivity
%  matricies in elemLst are adjusted for the new node numbering
%    
  
  numLst=length(elemLst);
  
  usedNodes=[];

  for n=1:numLst
    newNodes=unique(elemLst{n});
    if ( size(newNodes,2) ~= 1 )
      newNodes=newNodes';
    end
    usedNodes=[usedNodes;newNodes];
  end
  
  usedNodes=unique(usedNodes);
  numnode=size(usedNodes,1);
  
  nodeMap=zeros(size(node,1),1);
  nodeMap(usedNodes)=(1:numnode)';
   
  % remove unused nodes from node 
  node=node(usedNodes,:);  
  
  % renumber connectivity lists
  for n=1:numLst
    conn=elemLst{n};
    
    for e=1:size(conn,1)
      conn(e,:)=nodeMap(elemLst{n}(e,:))';
    end
    
    elemLst{n}=conn;
  end
  
  % Done !!!
  
    