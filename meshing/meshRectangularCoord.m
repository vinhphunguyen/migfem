function gcoord=meshRectangularCoord(L,D,numx,numy)
dx=L/numx;          % the length of side of 1 element in x-axis
dy=D/numy;          % the length of side of 1 element in y-axis
gcoord=[];
for j=1:numy+1
  for i=1:numx+1
      gcoord=[gcoord; (i-1)*dx (j-1)*dy;];
  end
end

