function conn=tricheck(node,conn,verbose)

% FUNCTION
%
%   conn=tricheck(node,conn,verbose)
%
% This function check wether a triangle has a negative Jacobian, and if
% so reorders it so that the the Jacobian is positive.

if ( nargin==2 )
  verbose=0;
end

if ( size(node,2)==3 )
  node=node(:,1:2);
end

count=0;

for e=1:size(conn,1)
  
  sctr=conn(e,:);
  [N,dNdxi]=lagrange_basis('T3',[1/3 1/3]);
  detJ=det(node(sctr,:)'*dNdxi);
  
  if ( detJ < 0 )
    %disp(['NEGATIVE JACOBIAN IN ELEMENT ',num2str(e)])
    conn(e,:)=fliplr(sctr);
    count=count+1;
  elseif ( detJ == 0 )
    disp(['ZERO JACOBIAN IN ELEMENT ',num2str(e),' CANNOT FIX'])
  end
end

if ( verbose )
  disp(['TRICHECK FOUND ',num2str(count),' NEGATIVE JACOBIANS, ALL FIXED'])
end