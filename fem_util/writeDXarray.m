function writeDXarray(A,fid,objectID,type)

% function writeDXarray(A,fid,objectID,type)
%
% This function writes data A in the openDx native format to the file given by
% its file id fid. objectId is the id given to the section and type is the data
% type ie. int, double, float..
%

rows = size(A,1);
cols = size(A,2);

if ( cols==1 )
  dataRank=0;
else
  dataRank=1;
end

dataItems=rows;
dataShape=cols;

if ( ischar(objectID) )
  objectID=['"',objectID,'"'];
else
  objectID=num2str(objectID);
end

% *************************** W R I T E   H E A D E R  **********************
% OBJECT
fprintf(fid,'object ');
fprintf(fid,objectID);

% CLASS
fprintf(fid,' class array ');

% TYPE
if ( strcmp(type,'double') )
  fprintf(fid,' type double ');
  dataFormat='%10.5b';
elseif ( strcmp(type,'int') )
  fprintf(fid,' type int ');
  dataFormat='%7d';
else
  fprintf(fid,' type float ');
  dataFormat='%10.5f';
end

% RANK
fprintf(fid,' rank ');
fprintf(fid,'%5d',dataRank);

% SHAPE
if ( dataRank~=0 )
    fprintf(fid,' shape ');
    fprintf(fid,'%3d',dataShape);
end

% ITEMS
fprintf(fid,' items ');
fprintf(fid,'%5d',dataItems);

% DATA LOCATION
fprintf(fid,' data follows ');
fprintf(fid,'\n');

% ****************************  W R I T E   D A T A  ***********************
for i=1:rows
    for j=1:cols
        
        fprintf(fid,dataFormat,A(i,j));
        fprintf(fid,' ');
        
    end
    fprintf(fid,'\n');
end
