function plot_field(X,connect,elem_type,field)
  
% function plot_field(X,connect,elem_type,field)
 
if ( nargin == 4 )
  nodesoff=0;
end
  
if ( size(field) == size(connect) )
  elementalField=1;
else
  elementalField=0;
end

% fill X if needed
if (size(X,2) < 3)
   for c=size(X,2)+1:3
      X(:,c)=[zeros(size(X,1),1)];
   end
end

holdState=ishold;
hold on

% plot elements
if     ( strcmp(elem_type,'Q9') )  % Q9 element
  ord=[1,5,2,6,3,7,4,8,1];
elseif ( strcmp(elem_type,'T3') )  % T3 element
  ord=[1,2,3,1];
elseif ( strcmp(elem_type,'T4') )  % T4 element
  ord=[1,2,3,1];
elseif ( strcmp(elem_type,'T6') )  % T6 element
  ord=[1,4,2,5,3,6,1];
elseif ( strcmp(elem_type,'Q4') )  % Q4 element
  ord=[1,2,3,4,1];
elseif ( strcmp(elem_type,'Q8') )  % Q8 element
   ord=[1,5,2,6,3,7,4,8,1];
elseif ( strcmp(elem_type,'L2') )  % L2 element
  ord=[1,2];   
elseif ( strcmp(elem_type,'L3') )  % L3 element
  ord=[1,3,2];   
end

for e=1:size(connect,1)
  
   xpt=X(connect(e,ord),1);
   ypt=X(connect(e,ord),2);      
   zpt=X(connect(e,ord),3);
   
   if ( elementalField )
     fpt=field(e,ord);
   else
     fpt=field(connect(e,ord));
   end
   
   fill3(xpt,ypt,zpt,fpt)
end

shading interp
axis equal
      
if ( ~holdState )
  hold off
end
