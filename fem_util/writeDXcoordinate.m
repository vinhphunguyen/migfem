function writeDXcoordinate(node,fid,objectID)

% function writeDXcoordinate(node,fid,objectID)
%
% Writes nodal coordinate matrix 

fprintf(fid,'# A NODE COORDINATE MATRIX\n');
writeDXarray(node,fid,objectID,'float')

fprintf(fid,'\n');