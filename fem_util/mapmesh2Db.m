function a=mapmesh2D(e1,e2,e3,e4)

% MAPMESH(e1,e2,e3,e4) Creates a mapped mesh along 2D surface
%
%  

if ( size(e1,1) ~= size(e3,1) |  size(e2,1) ~= size(e4,1) )
  disp('ERROR: opposite edges must be of same length')
end

if ( nargin == 3 )   % TRIANGULAR SURFACE

  

else                 % RULED SURFACE  
  
  nx=size(e1,1);
  ny=size(e2,1);

  % get straight edged reference grid
  p1=[e1(1,:)+e4(1,:)]/2;
  p2=[e1(nx,:)+e2(1,:)]/2;
  p3=[e2(ny,:)+e3(nx,:)]/2;
  p4=[e3(1,:)+e4(ny,:)]/2;
  
  strmesh=square_node_array(p1,p2,p3,p4,nx,ny);
  
  strx=zeros(nx,ny);
  stry=zeros(nx,ny);
  % put strmesh in ordered grid form
  for i=1:nx
    for j=1:ny
      strx(i,j)=strmesh((j-1)*nx+i,1);
      stry(i,j)=strmesh((j-1)*nx+i,2);
    end
  end
  clear strmesh;
  
  % get node push vectors along edges
  for i=1:nx
    dx(i,1)=e1(i,1)-strx(i,1);
    dy(i,1)=e1(i,2)-stry(i,1);
    dx(i,ny)=e3(i,1)-strx(i,ny);
    dy(i,ny)=e3(i,2)-stry(i,ny);
  end
  
  for j=1:ny
    dx(1,j)=e4(j,1)-strx(1,j);
    dy(1,j)=e4(j,2)-stry(1,j);
    dx(nx,j)=e2(j,1)-strx(nx,j);
    dy(nx,j)=e2(j,2)-stry(nx,j);
  end
  
  %interpolate nodal push vectors for interior nodes
  for i=2:nx-1
    for j=2:ny-1
      
      dx(i,j)=( (dx(1,j)*(nx-i)+dx(nx,j)*(i-1))/(nx-1)...
	       +(dx(i,1)*(ny-j)+dx(i,ny)*(j-1))/(ny-1) )/2;
      
      dy(i,j)=( (dy(1,j)*(nx-i)+dy(nx,j)*(i-1))/(nx-1)...
	       +(dy(i,1)*(ny-j)+dy(i,ny)*(j-1))/(ny-1) )/2;
      
    end
  end

  % push the straight mesh
  x=strx+dx;
  y=stry+dy;
  clear strx; clear stry; clear dx; clear dy;
  
  % unorder the mesh
  for i=1:nx
    for j=1:ny
      a((j-1)*nx+i,1)=x(i,j);
      a((j-1)*nx+i,2)=y(i,j);
    end
  end
  
  %clf
  %scatter(a(:,1),a(:,2))  

end
