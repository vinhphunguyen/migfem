function [varName,var,conn,time]=tecplotin(file)
%file='/home/jack/FreeSurface2D/step.plt';
% function [variable names,variable values,connectivity,time]=tecplotin(file)
%
% This function reads in a Tecplot data format.  It assumes that we are
% reading in FEPOINT data format.

% open the file
fid=fopen(file,'r');

if ( fid == -1 )  
  disp('ERROR: Was not able to open file')
else
  
  zoneNum=0;
  time=[];
  
  while 1
    
    line=fgetl(fid);            % read line
    if ~isstr(line), break, end % check if EOF
    
    if ( length(line)>=4 )
      if ( line(1)~='#' )  % then it is not a comment line
        
        testSeg=strtok(line,' =,'); 
        
        switch testSeg
          
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        case 'TITLE'                            % read in the TITLE line
          
          is=findstr(line,'"');
          
          if ( length(is)==2 )
            title=line(is(1)+1:is(2)-1);
          else
            disp('ERROR:syntax error in tecplot file')
          end
          
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        case 'VARIABLES'                      % read in the VARIABLE line
          
          is=findstr(line,'"');
          
          if ( mod(length(is),2) == 0 )
            
            numVar=length(is)/2;
            
            for i=1:numVar
              varName{i}=line(is(2*i-1)+1:is(2*i)-1);
            end
            
          else
            disp('ERROR:syntax error in tecplot file')
          end
          
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        case 'ZONE'                            % read in the ZONE info
          
          zoneNum=zoneNum+1;  
          line=line(5:length(line));
          
          while 1
            
            seg=strtok(line,',');
            [param,val]=strtok(seg,'=');
            param=deblank(strjust(param,'left'));
            val=val(2:length(val));
            
            if ( strcmp(param,'N') | strcmp(param,'E') )
              eval([seg,';']);
            elseif ( strcmp(param,'ET') )
              ET=val;
            elseif ( strcmp(param,'F') )
              F=val;
            elseif ( strcmp(param,'T') )
              is=findstr(val,'"');
              T=val(is(1)+1:is(2)-1);
              [a,b]=strtok(T,'=');
              t=str2num(b(2:length(b)));
            elseif ( strcmp(param,'D') )
              [a,b]=strtok(line,'(');
              a=strtok(b,')');
              [b,a]=strtok(a,'(');
              ndouble=0;
              while 1  % WORK on this
                ndouble=ndouble+1;
                
                break;
              end
            end
            
            % remove seg from line
            i=findstr(line,seg);
            start=findstr(line,seg)+length(seg)+1;
            stop=length(line);
            
            if ( start >= stop ) % test for eol
              break
            else
              line=line(start:stop);
            end
            
          end % of while loop
          time=[time;t]; 
          
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          % read zone node data
          val=zeros(N,numVar);
          n=0;
          while ( n < N )
            line=fgetl(fid);
            if ( length(line)~=0 )
              n=n+1;
              val(n,:)=str2num(line);
            end
          end
          
          var{zoneNum}=val;
          
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          % read connectivity data
          
          if ( zoneNum==1 )
            if ( strcmp(ET,'QUADRILATERAL') )
              r=4;
            else
              r=3;
            end
            
            conn=zeros(E,r);
            e=0;
            while ( e < E )
              line=fgetl(fid);
              if ( length(line)~=0 )
                e=e+1;
                conn(e,:)=str2num(line);
              end
            end
            
          end % of if first zoneNum
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          
        otherwise
          
          % do nothin
          
        end % of switch structure
        
      end % of if not a comment condition
    end  % or blank line
    
  end % of while loop
  
end  
