function [node,elemLst]=make_cross_mesh(pt1,pt2,numx,numy,ratiox,ratioy)

% function make_cross_mesh
%  
% [node,elemLst]=make_cross_mesh(p1,p2,numx,numy,ratiox,ratioy)
%
%    Makes a triangular cross mesh
%
 
  if ( nargin == 4 )
    ratiox=1;
    ratioy=1;
  end
  
  %pt1=[0,0]; pt2=[10,10]; numx=11; numy=11; ratiox=0.5; ratioy=0.5;
  pt1=pt1; pt3=pt2; pt2=[pt3(1) pt1(2)]; pt4=[pt1(1) pt3(2)];
  
  node=square_node_array(pt1,pt2,pt3,pt4,numx,numy,ratiox,ratioy);
  %clf; plot(node(:,1),node(:,2),'r*')
  %dx=(pt3(1)-pt1(1))/(numx-1);
  %dy=(pt3(2)-pt1(2))/(numy-1);
  
  %mpt1=pt1+[  dx/2,  dy/2 ]; 
  %mpt2=pt2+[ -dx/2,  dy/2 ];
  %mpt3=pt3+[ -dx/2,  -dy/2 ]; 
  %mpt4=pt4+[  dx/2,  -dy/2 ];
  
  %n=numx*numy;
  
  %mpt1=(node(1,:)+node(2,:)+node(numx+1,:)+node(numx+2,:))/4;
  %mpt2=(node(numx-1,:)+node(numx,:)+node(2*numx-1,:)+node(2*numx,:))/4;
  %mpt3=(node(n-1,:)+node(n,:)+node(n-numx,:)+node(n-numx-1,:))/4;
  %mpt4=(node(n-numx+1,:)+node(n-numx+2,:)+node(n-2*numx+1,:)+...
  %          node(n-2*numx+2,:))/4;
 
  %midnode=square_node_array(mpt1,mpt2,mpt3,mpt4,numx-1,numy-1,ratiox,ratioy);
  midnode=zeros((numx-1)*(numy-1),2);
  n=1;
  for j=1:(numy-1)
    for i=1:(numx-1)
      n1=numx*(j-1)+i;
      n2=n1+1;
      n4=numx*j+i;
      n3=n4+1;
      midnode(n,:)=mean(node([n1 n2 n3 n4],:));
      n=n+1;
    end
  end
  
  %hold on; plot(midnode(:,1),midnode(:,2),'b*')
  node=[node;midnode];
  
  N=numx*numy;
  
  connBot   = [ 1:(numx-1);
    2:numx]';
  connRight = [ numx:numx:numx*(numy-1);
    2*numx:numx:N]';
  connTop   = [ N:-1:(numx*(numy-1)+2);
    (N-1):-1:(numx*(numy-1)+1)]';
  connLeft  = [ (numx*numy-numx+1):-numx:numx;
                (numx*(numy-2)+1):-numx:1]';
  
  numCells = (numx-1)*(numy-1); 
  
  mdPt=N+1;
  connPtrn = [     1      2  mdPt;
	           2 numx+2  mdPt; 
              numx+2 numx+1  mdPt; 
	      numx+1      1  mdPt ];	
  
  conn=zeros(4*numCells,3);
  e=0;  
  for c=1:(numy-1)
    for r=1:(numx-1)
    
      conn( e*4 + [1 2 3 4], : ) = connPtrn;
      connPtrn=connPtrn+1;
      e=e+1;
      
    end
    connPtrn(:,1:2)=connPtrn(:,1:2)+1;   
  end
  
  elemLst{1}=connBot;
  elemLst{2}=connRight;
  elemLst{3}=connTop;
  elemLst{4}=connLeft;
  elemLst{5}=conn;