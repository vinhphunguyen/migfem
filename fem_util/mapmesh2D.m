function a=mapmesh2D(e1,e2,e3,e4)

% MAPMESH(e1,e2,e3,e4) Creates a mapped mesh along 2D surface
%
%  

if ( size(e1,1) ~= size(e3,1) |  size(e2,1) ~= size(e4,1) )
  disp('ERROR: opposite edges must be of same length')
end

nx=size(e1,1);
ny=size(e2,1);

if ( nargin == 3 )   % TRIANGULAR SURFACE

  

else                 % RULED SURFACE  
   
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % interpolate cubic serindipity mapping nodes from boundary data
  pt1=[e1(1,:)+e4(1,:)]/2;
  pt2=[e1(nx,:)+e2(1,:)]/2;
  pt3=[e2(ny,:)+e3(nx,:)]/2;
  pt4=[e3(1,:)+e4(ny,:)]/2;

  % on lower edge
  lseg=zeros(nx,1);
  for i=2:nx
    lseg(i)=lseg(i-1)+norm(e1(i,:)-e1(i-1,:));
  end

  L1=lseg(nx)/3; L2=2*lseg(nx)/3;
  temp=[interp1(lseg,e1(:,1),[L1,L2],'cubic');
        interp1(lseg,e1(:,2),[L1,L2],'cubic')]';   

  pt5=temp(1,:);  pt6=temp(2,:);
  E1=lseg*2/lseg(nx)-1;   % boundary nodes in parent domain

  % on right edge
  lseg=zeros(ny,1);
  for i=2:ny
    lseg(i)=lseg(i-1)+norm(e2(i,:)-e2(i-1,:));
  end

  L1=lseg(ny)/3; L2=2*lseg(ny)/3;
  temp=[interp1(lseg,e2(:,1),[L1,L2],'cubic');
        interp1(lseg,e2(:,2),[L1,L2],'cubic')]';   

  pt7=temp(1,:);  pt8=temp(2,:);
  E2=lseg*2/lseg(ny)-1;   % boundary nodes in parent domain

  % on top edge
  lseg=zeros(nx,1);
  for i=2:nx
    lseg(i)=lseg(i-1)+norm(e3(i,:)-e3(i-1,:));
  end

  L1=lseg(nx)/3; L2=2*lseg(nx)/3;
  temp=[interp1(lseg,e3(:,1),[L1,L2],'cubic');
        interp1(lseg,e3(:,2),[L1,L2],'cubic')]';   

  pt10=temp(1,:);  pt9=temp(2,:);
  E3=lseg*2/lseg(nx)-1;   % boundary nodes in parent domain    

  % on left edge
  lseg=zeros(ny,1);
  for i=2:ny
    lseg(i)=lseg(i-1)+norm(e4(i,:)-e4(i-1,:));
  end

  L1=lseg(ny)/3; L2=2*lseg(ny)/3;
  temp=[interp1(lseg,e4(:,1),[L1,L2],'cubic');
        interp1(lseg,e4(:,2),[L1,L2],'cubic')]';   

  pt12=temp(1,:);  pt11=temp(2,:);
  E4=lseg*2/lseg(ny)-1;   % boundary nodes in parent domain 

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % setup interior node grid on parent domain
  x=zeros(nx,ny); y=zeros(nx,ny);

  X=[pt1(1); pt2(1); pt3(1); pt4(1); pt5(1); pt6(1);
      pt7(1); pt8(1); pt9(1); pt10(1); pt11(1); pt12(1)]; 
 
  Y=[pt1(2); pt2(2); pt3(2); pt4(2); pt5(2); pt6(2);
      pt7(2); pt8(2); pt9(2); pt10(2); pt11(2); pt12(2)]; 


  for i=2:nx-1
    for j=2:ny-1

      A=[E4(j)-E2(j) 2;2 E1(i)-E3(i)];
      b=[E2(j)+E4(j);E1(i)+E3(i)];
      inter=(A\b);
      xi=inter(1); eta=inter(2);
   
      N(1)=(1-xi)*(1-eta)*(9*(xi^2+eta^2)-10)/32;     % mapping functions
      N(2)=(1+xi)*(1-eta)*(9*(xi^2+eta^2)-10)/32; 
      N(3)=(1+xi)*(1+eta)*(9*(xi^2+eta^2)-10)/32;  
      N(4)=(1-xi)*(1+eta)*(9*(xi^2+eta^2)-10)/32;

      N(5)=9/32*(1-eta)*(1-xi^2)*(1-3*xi);
      N(6)=9/32*(1-eta)*(1-xi^2)*(1+3*xi);
       
      N(7)=9/32*(1+xi)*(1-eta^2)*(1-3*eta);
      N(8)=9/32*(1+xi)*(1-eta^2)*(1+3*eta);
  
      N(9)= 9/32*(1+eta)*(1-xi^2)*(1+3*xi);
      N(10)=9/32*(1+eta)*(1-xi^2)*(1-3*xi); 
       
      N(11)=9/32*(1-xi)*(1-eta^2)*(1+3*eta);
      N(12)=9/32*(1-xi)*(1-eta^2)*(1-3*eta);        

      x(i,j)=N*X; y(i,j)=N*Y;
      
    end
  end

  x(:,1)  = e1(:,1);
  x(nx,:) = e2(:,1)';
  x(:,ny) = e3(:,1);
  x(1,:)  = e4(:,1)';

  y(:,1)  = e1(:,2);
  y(nx,:) = e2(:,2)';
  y(:,ny) = e3(:,2);
  y(1,:)  = e4(:,2)';
   
  a=zeros(nx*ny,2);
  for j=1:ny
    for i=1:nx
      a((j-1)*nx+i,1)=x(i,j); 
      a((j-1)*nx+i,2)=y(i,j);
    end
  end 

end
