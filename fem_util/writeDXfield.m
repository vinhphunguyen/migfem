function writeDXfield(fid,fieldID,coordID,connID,dataID)

% function writeDXfield(fid,fieldID,coordID,connID,dataID)
%
% Writes a FEM field
%
% example:
% node=[0 0;1 0;2 0;0 1;1 1;2 1;0 2;1 2;2 2];
% conn=[1     2     4;
%       2     5     4;
%      2     6     5;
%      2     3     6;
%      4     8     7;
%      4     5     8;
%      5     6     8;
%      6     9     8];
%  d1=[0 0 0 1 1 1 2 2 2]';
%  d2=d1+1;
%  d3=d2+0.5;
%  
%  filename='datafile.dx';
%  fid=fopen(filename,'w');
%  writeDXcoordinate(node,fid,1);
%  writeDXconnectivity(conn,fid,2,'T3');
%  
%  writeDXdata(d1,fid,3);
%  writeDXfield(fid,103,1,2,3);
%  
%  writeDXdata(d2,fid,4);
%  writeDXfield(fid,104,1,2,4);
%  
%  writeDXdata(d3,fid,5);
%  writeDXfield(fid,105,1,2,5);

if ( ischar(fieldID) )
  fieldID=['"',fieldID,'"'];
else
  fieldID=num2str(fieldID);
end

if ( ischar(coordID) )
  coordID=['"',coordID,'"'];
else
  coordID=num2str(coordID);
end

if ( ischar(connID) )
  connID=['"',connID,'"'];
else
  connID=num2str(connID);
end

if ( ischar(dataID) )
  dataID=['"',dataID,'"'];
else
  dataID=num2str(dataID);
end

fprintf(fid,'# FIELD DEFINITION\n');
fprintf(fid,'object ');
fprintf(fid,fieldID);
fprintf(fid,' class field\n');

fprintf(fid,'component \"positions\" value ');
fprintf(fid,coordID);

fprintf(fid,'\ncomponent \"connections\" value ');
fprintf(fid,connID);

fprintf(fid,'\ncomponent \"data\" value ');
fprintf(fid,dataID);

fprintf(fid,'\n\n');

