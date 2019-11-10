function [splitNode,splitConn]=splitElement(node,conn,elemType,nnn)

switch (elemType)
  
case 'T3'
  nnN=[0 0 1;1 0 0;0 1 0;0.5 0 0.5;.5 .5 0;0 .5 0.5];
  nnN=nnN(:,[2 3 1]);
  splitNode=nnN*node;
  
  n1=conn(1); n2=conn(2); n3=conn(3);
  splitConn=[n1 nnn nnn+2;nnn n2 nnn+1;nnn+2 nnn+1 n3;nnn nnn+1 nnn+2];
  
otherwise
  disp(['ERROR: ',elemType,' not defined in splitElement function'])  
end