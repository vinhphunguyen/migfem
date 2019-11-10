function element=make_elem(node_pattern,num_u,num_v,inc_u,inc_v)

% function element=make_elem(node_pattern,num_u,num_v,inc_u,inc_v)
%
% creates a connectivity list

if ( nargin < 5 )
   disp(['Not enough parameters specified for make_elem function'])
end

inc=[zeros(1,size(node_pattern,2))];
e=1;
element=zeros(num_u*num_v,size(node_pattern,2));

for row=1:num_v
   for col=1:num_u
      element(e,:)=node_pattern+inc;
      inc=inc+inc_u;
      e=e+1;
   end
   inc=row*inc_v;
end
