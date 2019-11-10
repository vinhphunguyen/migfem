function a=linemesh(p1,p2,n,ratio)

% LINEMESH(P1,P2,N,R) Creates a mesh along a line
%
%  

if ( nargin < 4 )
  ratio=1;
end

ratio=abs(ratio);

if ( ratio==1 )
  a=[linspace(p1(1),p2(1),n);linspace(p1(2),p2(2),n)]';
else

  rv=ratio^(1/(n-2));   % get scaling
  w(1)=0;
  d=1;
  for i=2:n
    w(i)=w(i-1)+d;
    d=d/rv;
  end
  w=w'/w(n);  

  a=(ones(n,1)-w)*p1+w*p2;  

end
