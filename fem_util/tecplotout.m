function tecplotout(element,elemtype,node,varname,value,filename,append,time)
  
% function tecplotout(element,elemtype,node,varname,value,filename)
%
% This function writes out a tecplot FEPOINT style data file  

if ( nargin==6 )
  append=0;
  time=0.0;
end
  
if ( nargin==7 )
  time=0.0;
end
  
if ( append )
  fid=fopen(filename,'a');
else
  fid=fopen(filename,'w');
end

if ( fid == -1 )
  disp('ERROR:Could not open file')
end

numZones=length(value);
numVar=length(varname);
numNodes=length(node);
numElem=length(element);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% write file header info

% ****** Write the VARIABLES section
varStr=['VARIABLES = "',varname{1},'"'];
for v=2:numVar
  varStr=[varStr,', "',varname{v},'"'];
end

fprintf(fid,[varStr,'\n\n']);


% ****** Write the ZONE(s) section
switch elemtype
 case 'T3'
  etStr='TRIANGLE';
 case 'Q4'
  etStr='QUADRILATERAL';
 otherwise
  etStr='TRIANGLE';
end

fprintf(fid,['ZONE N=',num2str(numNodes),', E=',num2str(numElem),...
      ', F=FEPOINT, ET=',etStr,', T="t=',num2str(time),'"\n\n']);

% get format strings
valFormat=' %8.4f';
for i=2:numVar
  valFormat=[valFormat,' %8.4f'];
end
valFormat=[valFormat,' \n'];

connFormat=' %5i';
for i=2:size(element,2);
  connFormat=[connFormat,' %5i'];
end
connFormat=[connFormat,' \n'];

if ( iscell(value) )  
  
  fprintf(fid,valFormat,(value{1})'); % print values for ZONE z
  fprintf(fid,'\n');
  fprintf(fid,connFormat,element');        % print connectivity
  fprintf(fid,'\n\n');
  
  for z=2:numZones
  
    fprintf(fid,['ZONE N=',num2str(numNodes),', E=',num2str(numElem),...
        ', F=FEPOINT, ET=',etStr,' D=(FECONNECT), T="t=',num2str(time),'"\n\n']);
  
    fprintf(fid,valFormat,(value{z})');   % print values for ZONE z
    fprintf(fid,'\n\n');
    
  end  % of ZONE loop
  
else
  
  fprintf(fid,valFormat,value');   % print values for ZONE z
  fprintf(fid,connFormat,element');    % print connectivity
  
end

% close file
status=fclose(fid);
