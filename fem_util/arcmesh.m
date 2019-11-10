function a=arcmesh(p1,p2,r,n,ratio)

% ARCMESH(P1,P2,N,R) Creates a mesh along an arc
%
%  

if ( nargin < 5 )
  ratio=1;
end

ratio=abs(ratio);

d=norm(p2-p1);
l=sqrt(r^2-d^2/4);

if ( r > 0 )
  center=(p1+p2)/2-l*[p2(2)-p1(2),p1(1)-p2(1)]/d;
else
  center=(p1+p2)/2+l*[p2(2)-p1(2),p1(1)-p2(1)]/d;
end
  
r1=p1-center;
r2=p2-center;

theta1=atan2(r1(2),r1(1));
theta2=atan2(r2(2),r2(1));

if ( ratio==1 )

  theta=(linspace(theta1,theta2,n))';
  
else

  rv=ratio^(1/(n-2));   % get scaling
  w(1)=0;
  d=1;
  for i=2:n
    w(i)=w(i-1)+d;
    d=d/rv;
  end
  w=w'/w(n);  

  theta=(ones(n,1)-w)*theta1+w*theta2;  

end

a=ones(n,1)*center+r*[cos(theta),sin(theta)];
a(1,:)=p1;
a(n,:)=p2;
