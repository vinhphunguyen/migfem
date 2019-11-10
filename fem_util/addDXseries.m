function addDxseries(filename,seriesID,fieldIDs,timeIDs)

% function addDxseries(filename,seriesname,fieldIDs,timeIDs)
%
% adds a sereis to a dx data file

if ( nargin < 4 )  
  timeIDs=0:(length(fieldIDs)-1);
end

fid=fopen(filename,'a');

writeDXseries(fid,seriesID,fieldIDs,timeIDs);