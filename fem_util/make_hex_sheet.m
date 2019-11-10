numx=10;
numy=10;

nodes=zeros((2*numx+1)*(numy+1),2);
n=1;

for j=1:numy
  for i=1:numx
    
    node(n,:)=(
    
    n=n+1;
  end
end



cellCenter=square_node_array([0,0],[numx-1,0],[numx-1,numy-1],[0,numy-1]...
			     ,numx,numy);
cellCenterXY=cellCenter;
cellCenterXY(:,1)=cellCenter(:,1)+(1/2)*cellCenter(:,2);
cellCenterXY(:,2)=(sqrt(3)/2)*cellCenter(:,2);
figure(1); clf; plot(cellCenterXY(:,1),cellCenterXY(:,2),'bd')

node=zeros(numx*numy*3,2);
conn=zeros(numx*numy*3,2);

for c=1:numy*numy
  node( c*3-[2 1 0], : )=[ cellCenter(c,:)+[ -1/3, -1/3 ];
		           cellCenter(c,:)+[  1/3, -2/3 ];
		           cellCenter(c,:)+[  2/3, -1/3 ] ];
  cup=c+numx;
  conn( c*3-[2 1 0], : )=[ 3*c-2 3*c-1;
		           3*c-1   3*c;
		           3*c   3*cup-1 ];
   
end

nodeXY=node;
nodeXY(:,1)=node(:,1)+(1/2)*node(:,2);
nodeXY(:,2)=(sqrt(3)/2)*node(:,2);
figure(1); hold on; plot_mesh(nodeXY,conn,'R2');
