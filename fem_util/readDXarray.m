function A=readDXarray(filename,id)

% function data=readDXarray(filename,id)
%
% reads a data field from an open dx file

% open file
% clear
% fid=fopen('surf_tension_test.dx','r');
% id = 20000;

if ( ischar(id) )
  id=['"',id,'"'];
else
  id=num2str(id);
end

fid = fopen(filename,'r');
if ( fid == -1 )
  disp(['ERROR: COULD NOT OPEN "',filename,'"'])
  return
end

% fid section
while 1
  
  line=fgetl(fid);            % read line
  if ~isstr(line), break, end % check if EOF
 
  if ( length(line) == 0 )
    line ='#';
  end
    
  if ( line(1)~='#' )  % then it is not a comment line
    
    [testSeg,restOfLine] = strtok(line,' ,'); 
    
    if strcmp(testSeg,'object')
    
      [oid,restOfLine] = strtok(restOfLine,' ,'); 
      
      if ( ~ischar(oid) )
        oid=num2str(oid);
      end
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if ( strcmp(oid,id) )  % found it !!!
        
        [testSeg,restOfLine] = strtok(restOfLine,' ,'); 
        % find rank 
        while ( ~strcmp(testSeg,'rank') )
          [testSeg,restOfLine] = strtok(restOfLine,' ,'); 
        end
        [rank,restOfLine] = strtok(restOfLine,' ,'); 
        rank=str2num(rank);
        
        if ( rank~=0 ) % find shape
          
          [testSeg,restOfLine] = strtok(restOfLine,' ,'); 
          while ( ~strcmp(testSeg,'shape') )
            [testSeg,restOfLine] = strtok(restOfLine,' ,'); 
          end
          [shape,restOfLine] = strtok(restOfLine,' ,'); 
          shape=str2num(shape);
          
        else
          shape = 1;
          rank  = 1;
        end
        
        % find items
        [testSeg,restOfLine] = strtok(restOfLine,' ,'); 
        while ( ~strcmp(testSeg,'items') )
          [testSeg,restOfLine] = strtok(restOfLine,' ,'); 
        end
        [items,restOfLine] = strtok(restOfLine,' ,'); 
        items=str2num(items);
        
        % read array
        A=zeros( items, rank*shape );
        
        n=0;
        while ( n < items )
          line=fgetl(fid);
          if ( length(line)~=0 )
            n=n+1;
            A(n,:)=str2num(line);
          end
        end
        
        break
      end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
    end
    
  end % if comment line
  
end
fclose(fid)