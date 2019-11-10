function [tri,X,Xi]=getSubElements2(node,ne)

%   function [tri,X,Xi]=getSubElements2(node,ne)
%
  
  
% node=[0 0;1 0;0 1];
% d=[1 2 3 2 1 5]';
% phi=[-1 2 .5]';
% ne=25;

if ( nargin == 3 )
  ne=25;
end

% subdivide elements
% get nodes and nodal values
numnode=sum(ne:-1:1);
X=zeros(numnode,2);
Xi=zeros(numnode,2);
del=1/(ne-1);
n=1;
eta=0.0;
for j=1:ne   
  xi=0.0;
  for i=1:(ne-j+1)
    N=lagrange_basis('T3',[xi,eta]);
    X(n,:)=N'*node;
    Xi(n,:)=[xi,eta];
    
    n=n+1;
    xi=xi+del;
  
  end
  eta=eta+del;
  
end

% get elements
tri=delaunay(X(:,1),X(:,2));
