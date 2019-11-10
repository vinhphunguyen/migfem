function writeDXconnectivity(conn,fid,objectID,elemType)

% function writeDXconnectivity(conn,fid,objectID,elemType)
%
% Writes element connectivity matrix 

fprintf(fid,'# AN ELEMENT CONNECTIVITY MATRIX\n');
writeDXarray(conn-1,fid,objectID,'int');

if ( strcmp(elemType,'L2') )
  fprintf(fid,'attribute \"element type\" string \"lines\" \n');
elseif ( strcmp(elemType,'T3') )
  fprintf(fid,'attribute \"element type\" string \"triangles\" \n');
elseif ( strcmp(elemType,'Q4') )
  fprintf(fid,'attribute \"element type\" string \"quads\" \n');
elseif ( strcmp(elemType,'H4') )
  fprintf(fid,'attribute \"element type\" string \"tetrahedra\" \n');
elseif ( strcmp(elemType,'B8') )
  fprintf(fid,'attribute \"element type\" string \"cubes\" \n');
end

fprintf(fid,'attribute \"ref\" string \"positions\" \n\n');
