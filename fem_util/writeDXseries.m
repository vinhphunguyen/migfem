function writeDXseries(fid,objectID,fieldIDs,positions)

% function writeDXseries(fid,objectID,fieldIDs,positions)
%
% Writes a series
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
%  
%  writeDXseries(fid,1000,103:105,0:0.2:0.4);

N=length(fieldIDs);

if ( nargin==3 )
    positions=0:N-1;
end

if ( ischar(objectID) )
  objectID=['"',objectID,'"'];
else
  objectID=num2str(objectID);
end

fprintf(fid,'# SERIES DEFINITION\n');
fprintf(fid,'object ');
fprintf(fid,objectID);
fprintf(fid,' class series\n');

for n=1:N
    mem=n-1;
    val=fieldIDs(n);
    pos=positions(n);
    
    fprintf(fid,'member ');
    fprintf(fid,'%5d',mem);
    fprintf(fid,' value ');
    fprintf(fid,'%5d',val);
    fprintf(fid,' position ');
    fprintf(fid,'%8.3f',pos);
    fprintf(fid,'\n');
    
end