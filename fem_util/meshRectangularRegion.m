% =======================================================================
% function [node,element] = meshRectangularRegion(...
%    pt1, pt2, pt3, pt4, nnx,nny,elemType)
%
% elemType : element type 'Q4', 'Q9', 'PARTICLE'
%
% Generates a quadratleral array of nodes between the counterclockwise 
% ordering of nodes pt1 - pt4.  There are numnod_u nodes in the u direction
% (pt1 - pt2) and numnode_v nodes in the v direction (pt2 - pt3).  The 
% parameter uratio and vratio determint the nodal spacing along the u and v
% lines.  If no values of uratio and/or vratio are given values of unity are
% assumed which resulets in uniformed node spacng along the u an v directions.
% uratio and v ratio  are the ratio of the first node spacing to the last
% node spacing along th u of v direction respectivly (the first spacing
% occurs near pt 1 and the last near pt3. 
%
% GOAL :
% Generate the nodal coordinates (node) and connectivity (element) for a
% rectangular region given the four corners pti, the number of nodes in
% each direction (nnx,nny), the element type 
%
% =======================================================================
% Stephane Bordas 2004-November-30=2004.11.30
% LSC EPFL
% =======================================================================
function [node,element] = meshRectangularRegion(pt1, pt2, pt3, pt4, nnx,nny,elemType)

% obtain number of elements from number of nodes
numx = nnx-1;
numy = nny-1;

switch elemType
    case 'PARTICLE' % generate the mesh for particles
        % NOTE : for now, does just like a Q4
        % in actuality, the element returned by make_elem is just dummy and
        % is to be used nowhere. only kept so that element is returned and
        % all output arguments are defined when function returns.
        % SB.2004.12.01
        node=square_node_array(pt1,pt2,pt3,pt4,nnx,nny);
        inc_u=1;
        inc_v=nnx;
        node_pattern=[ 1 2 nnx+2 nnx+1 ];
        [element]=make_elem(node_pattern,numx,numy,inc_u,inc_v);
        
        
    case 'Q4'           % here we generate the mesh of Q4 elements
        node=square_node_array(pt1,pt2,pt3,pt4,nnx,nny);
        inc_u=1;
        inc_v=nnx;
        node_pattern=[ 1 2 nnx+2 nnx+1 ];
        [element]=make_elem(node_pattern,numx,numy,inc_u,inc_v);

    case 'Q9'           % here we generate a mesh of Q9 elements
        node=square_node_array(pt1,pt2,pt3,pt4,nnx,nny);
        inc_u=2;
        inc_v=2*nnx;
        node_pattern=[ 1 3 2*nnx+3 2*nnx+1 2 nnx+3 2*nnx+2 nnx+1 nnx+2 ];

        [element]=make_elem(node_pattern,numx,numy,inc_u,inc_v);
    case 'T3' % and last but not least T3 elements
        
        node=square_node_array(pt1,pt2,pt3,pt4,nnx,nny);
        node_pattern1=[ 1 2 nnx+1 ];
        node_pattern2=[ 2 nnx+2 nnx+1 ];
        inc_u=1;
        inc_v=nnx;
        element=[make_elem(node_pattern,numx,numy,inc_u,inc_v);
            make_elem(node_pattern,numx,numy,inc_u,inc_v)];
    otherwise
        error('For now, only PARTICLE, Q4, Q9 and T3 are supported by the mesh generator');
end



% END OF FUNCTION
