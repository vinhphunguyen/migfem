function conn=brickcheck(node,conn,verbose)

% FUNCTION
%
%   conn=brickcheck(node,conn,verbose)
%
% This function check wether a brick has a negative Jacobian, and if
% so reorders it so that the the Jacobian is positive.

if ( nargin==2 )
  verbose=0;
end

count=0;

for e=1:size(conn,1)
  
  sctr=conn(e,:);
  [N,dNdxi]=lagrange_basis('B8',[0 0 0]);
  detJ=det(node(sctr,:)'*dNdxi);
  
  if ( detJ < 0 )
    %disp(['NEGATIVE JACOBIAN IN ELEMENT ',num2str(e)])
    conn(e,:)=sctr([4 3 2 1 8 7 6 5]);
    count=count+1;
  elseif ( detJ == 0 )
    disp(['ZERO JACOBIAN IN ELEMENT ',num2str(e),' CANNOT FIX'])
  end
end

if ( verbose )
  disp(['BRICKCHECK FOUND ',num2str(count),' NEGATIVE JACOBIANS, ALL FIXED'])
end