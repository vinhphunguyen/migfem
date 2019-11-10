function [node,element,elemType]=msh2mlab(meshFile)

% function [node,connectivities,elemType]=msh2mlab('meshFile')  
%
% This file reads in a mesh file from gmsh and returns the appropirate input
% data structures for a matlab finite element code. 
%
%       'meshFile' - is the file name of the *.msh file ( includeing *.msh )
%
%       node - is an nodal coordinate matrix ( Nx3, where N is the number of nodes )
%
%       connectivities - is a cell array of connetivitiy matricies for each
%                        zone defined in the *.msh ( or *.geo ) file
%
%       elemType - is a vector of strings that tell what type of element is in 
%                  each zone ( we assume that only one type of element is 
%                  defined for each zone )
%
% Example:
%       plate.msh is a mesh file with three zones. Zone 1 is the mesh of
%       the domain interior, zone 2 is a mesh of the traction boundary
%       and zone 6 is a mesh of the displacement boundary the code to 
%       process this file would be 
%
%           [node,conns,elemType]=msh2mlab('plate.msh'); 
%           domain=conns{1};
%           tracBndy=conns{2};
%           dispBndy=conns{6};
%           dispBndyNodes=unique(dispBndy);
%

% open the file
%meshPath=''; %'/home/jack/Meshes/';
fid=fopen(meshFile,'r');
% define mesh data structures
pts=[]; seg={}; zon={}; node=[]; element={};

%* README: The 'msh' file format is the native output file format for
%  Gmsh. The file is divided in several sections (enclosed in $KEY and
%  $ENDKEY pairs). Two fields are important: $NOD/$ENDNOD defines the
%  nodes and $ELM/$ENDELM defines the elements.
%
%  The syntax is as follows:
%
%  $NOD
%  number-of-nodes
%  node-number x-coord y-coord z-coord 
%  ...
%  $ENDNOD
%
%  $ELM
%  number-of-elements
%  elm-number elm-type elm-region unused number-of-nodes node-numbers
%  ...
%  $ENDELM
%
%  All the syntactic variables stand for integers except x-coord,
%  y-coord and z-coord which stand for floating point values.  The
%  elm-type value defines the geometrical type for the element:
%  
%  elm-type: 
%  
%  1 Line (2 nodes, 1 edge). 
%  2 Triangle (3 nodes, 3 edges). 
%  3 Quadrangle (4 nodes, 4 edges). 
%  4 Tetrahedron (4 nodes, 6 edges, 4 facets). 
%  5 Hexahedron (8 nodes, 12 edges, 6 facets). 
%  6 Prism (6 nodes, 9 edges, 5 facets). 
%  7 Pyramid (5 nodes, 8 edges, 5 facets). 
%  15 Point (1 node). 
%
%  The elm-region value is the number of the physical entity to which
%  the element belongs. 
zoneID=0;
% read sections
while 1
  
  line=fgetl(fid);            % read line
  if ~isstr(line), break, end % check if EOF
  
  switch line                 % find seciton
    
  case '$PTS'               
    n=str2num(fgetl(fid));
    pts=zeros(n,6);
    
    for i=1:n
      pts(i,:)=str2num(fgetl(fid));
    end
    
  case '$SEG'              
    n=str2num(fgetl(fid));
    
    for i=1:n
      seg{i}=str2num(fgetl(fid));
    end           
    
  case '$ZON'             
    n=str2num(fgetl(fid));
    
    for i=1:n
      zon{i}=str2num(fgetl(fid));
    end                    
    
    zoneID=ones(length(zon),1);
    for i=1:length(zon);
      zoneID(i)=zon{i}(1);
      element{zoneID(i)}=[];
      elemType{i}='??';
    end
    
  case '$NOD'           
    numnode=str2num(fgetl(fid));
    node=zeros(numnode,3);
    
    for i=1:numnode
      nodeline=str2num(fgetl(fid));
      node(nodeline(1),:)=nodeline(2:4);
    end
    
  case '$NOE'           
    numnode=str2num(fgetl(fid));
    node=zeros(numnode,3);
    
    for i=1:numnode
      nodeline=str2num(fgetl(fid));
      node(nodeline(1),:)=nodeline(2:4);
    end
    
  case '$ELM'          
    n=str2num(fgetl(fid));
    
    for i=1:n
      temp=str2num(fgetl(fid));  % get element
      inZone=temp(3);            % find zone element is in
      z=find(ismember(zoneID,inZone));
      
      if ( isempty(z) ) % zone is not yet defined
        zoneID=[zoneID;inZone];
        z=size(zoneID,1);
        elemType{zoneID(z)}='??';
        element{zoneID(z)}=[];    
      end
      
      if isempty(element{z}) % first element inZone so set type
        switch temp(2)
        case 1
          elemType{zoneID(z)}='L2';
        case 2
          elemType{zoneID(z)}='T3';
        case 3
          elemType{zoneID(z)}='Q4';
        case 4
          elemType{zoneID(z)}='H4';
        case 5
          elemType{zoneID(z)}='B8';
        case 6
          elemType{zoneID(z)}='P6';
        case 15
          elemType{zoneID(z)}='P1';
        otherwise
          elemType{zoneID(z)}='??';
        end
      end
      
      % add element to that zone
      element{zoneID(z)}=[element{zoneID(z)};temp(6:5+temp(5))];
      
    end
    
  otherwise
    % skip line
  end
  
end

fclose(fid);
