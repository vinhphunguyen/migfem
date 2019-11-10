function [node,element,elemType] = sysmesh2mlab(meshFile)

% function [node,element,elemType] = sysmesh2mlab('meshFile')
%
% This file reads in a mesh file from SYSMESH (ESI) and returns the appropirate input
% data structures for a matlab finite element code.
%
%       'meshFile' - is the file name of the *.asc file
%
%       node - is an nodal coordinate matrix ( Nx3, where N is the number of nodes )
%
%       element - is a cell array of connetivitiy matricies for each
%                        zone defined in the *.msh ( or *.geo ) file
%
%       elemType - is a vector of strings that tell what type of element is in
%                  each zone ( we assume that only one type of element is
%                  defined for each zone )
%
%
%
%  All the syntactic variables stand for integers except x-coord,
%  y-coord and z-coord which stand for floating point values.  The
%  elm-type value defines the geometrical type for the element:
%
%  elm-type:
%
%  2003 Triangle (3 nodes, 3 edges).
%  2006 Triangle (6 nodes, 3 edges).
%  2004 Quadrangle (4 nodes, 4 edges).
%  2008 Quadrangle (8 nodes, 4 edges).
%  3004 Tetrahedron (4 nodes, 6 edges, 4 facets).
%  ...
%
%  The elm-region value is the number of the physical entity to which
%  the element belongs.

% Written by Vinh Phu NGUYEN, LTDS, ENISE, 25 May 2006

% open the file
fid = fopen(meshFile,'r');

% define mesh data structures
node=[]; element=[];
zon={};
zoneID=0;

% read sections
while 1

    line = fgetl(fid);            % read line
    if ~isstr(line), break, end % check if EOF

    switch line                 % find section
        % Node section
        case 'NODE'
            numnode = str2num(fgetl(fid));
            node = zeros(numnode,3);
            for i = 1:numnode
                nodeline = str2num(fgetl(fid));
                node(nodeline(1),:) = nodeline(7:9);
            end

        % Element section
        case 'ELEM'
            numelem = str2num(fgetl(fid));        % number of elements
            for i = 1:numelem
                temp = str2num(fgetl(fid));% get element
                switch temp(2)
                    case 2003
                        elemType  =  'T3';
                        nne       =   3 ;    % number of node per element
                    case 2006
                        elemType  =  'T6';
                        nne       =   6 ;    % number of node per element
                    case 2004
                        elemType  =  'Q4';
                        nne       =   4 ;    % number of node per element
                    case 2008
                        elemType  =  'Q8';
                        nne       =   8 ;    % number of node per element
                    otherwise
                        elemType  =  '??';
                end

                % add element to that zone
                element  =  [element;temp(6:5 + nne)];
            end

        otherwise
            % skip line
    end

end

fclose(fid);
