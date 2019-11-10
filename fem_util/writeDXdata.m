function writeDXdata(d,fid,objectID)

% function writeDXdata(d,fid,objectID)
%
% Writes nodal data

fprintf(fid,'# NODAL DATA ARRAY\n');
writeDXarray(d,fid,objectID,'float');

fprintf(fid,'attribute \"dep\" string \"positions\" \n\n');
