function [node,element,b1,b2,b3,b4,b5,b6] = ... 
            holeinplate(h,w,r,numr,numh,numw,ratio,type)
% function [node,element,b1,b2,b3,b4,b5,b6] = 
%                       holeinplate(h,w,r,numr,numh,numw,ratio,type)
%
% This function make a mesh of a half rectangle with a hole in it.  The
% input variables are as follows:
%   h: is the height of the rectangle    
%   w: is the width of the rectangle 
%   numr: is the number of nodes in the radial direction
%   numh: is the number of nodes in the half height direction
%   numw: is the number of nodes in the width direction
%   ratio: is a element edge ratio between an element on the hole to an
%          element on the outer edge
%   type: is the element type presently only 'Q4' and 'T3' are valid.
%
% The output variables are as follows:
%   node: is a node coordinate matrix
%   element: is an element connectivity matrix
%   b1-b6: are arrays of the boundary nodes on the edges 1-6 as numbered
%   below.  The nodes are numbered in counter-clockwise fashion
%
%                              b3 
%          --------------------------------------------
%          |                                          |
%          |                                          |
%          |                  ______                  |
%       b4 |                 /      \                 | b2
%          |                /   b5    \               |
%          |              /            \              |
%          |              |            |              |
%          ----------------            ----------------
%               b5                           b1 

% aspect angle
theta=atan2(h,w);

% get points
p1=[r,0];
p2=[w/2,0];
p3=[w/2,h/2];
p4=[cos(theta),sin(theta)]*r;
p5=[0,h/2];
p6=[0,r];

numu=numr;
numv=numh;
numv2=floor(numw/2)+1;

ratioh=(.5*(pi/2-theta)+1*theta)/pi*2;
ratiow=(1.5*(pi/2-theta)+1*theta)/pi*2;

% mesh boundaries
e1=linemesh(p1,p2,numu,ratio);
e2=linemesh(p2,p3,numv,ratioh);
e3=linemesh(p4,p3,numu,ratio);
e4=arcmesh(p1,p4,r,numv);

e5=linemesh(p3,p5,numv2,ratiow);
e6=linemesh(p6,p5,numu,ratio);
e7=arcmesh(p4,p6,r,numv2);

% make node mesh
node=mapmesh2D(e1,e2,e3,e4);
node=[node(1:numu*(numv-1),:);mapmesh2D(e3,e5,e6,e7)];
numnode=rows(node);

% mirror node mesh to left hand side
node_lhs=[-node(:,1),node(:,2)];
node=[node(1:numnode-numu,:);node_lhs];
numnode=rows(node);

% make elements
if ( type == 'Q4' )  %  Q4 ELEMENTS
  
  element1=make_elem([1 2 numu+2 numu+1],numu-1,numv+numv2-3,1,numu);
  offset=(numu)*(numv+numv2-2);
  element2=make_elem(offset+[1 2 numu+2 numu+1],numu-1,numv+numv2-2,1, ...
		   numu);
  n1=numu*(numv+numv2-3)+1;
  n2=numnode-numu+1;
  element3=make_elem([n1 n2 n2+1 n1+1],numu-1,1,1,numu);
  element=[element1;element2;element3];
  numelem=rows(element);

elseif ( type == 'T3' ) % T3 ELEMENT
  
  element1a=[make_elem([1 2 numu+1],numu-1,numv-1,1,numu);
	     make_elem([2 numu+2 numu+1],numu-1,numv-1,1,numu)];
  n0=numu*(numv-1);
  element1b=[make_elem(n0+[1 2 numu+2],numu-1,numv2-2,1,numu);
             make_elem(n0+[1 numu+2 numu+1],numu-1,numv2-2,1,numu)];
  element1=[element1a;element1b];
  offset=(numu)*(numv+numv2-2);
  element2a=[make_elem(offset+[1 2 numu+1],numu-1,numv+numv2-2,1,numu);
	  make_elem(offset+[2 numu+2 numu+1],numu-1,numv+numv2-2,1, ...
		    numu)];
  element2a=[make_elem(offset+[1 2 numu+1],numu-1,numv-1,1,numu);
	     make_elem(offset+[2 numu+2 numu+1],numu-1,numv-1,1,numu)];
  n0=numu*(numv-1)+offset;
  element2b=[make_elem(n0+[1 2 numu+2],numu-1,numv2-1,1,numu);
             make_elem(n0+[1 numu+2 numu+1],numu-1,numv2-1,1,numu)];
  element2=[element2a;element2b];
  
  n1=numu*(numv+numv2-3)+1;
  n2=numnode-numu+1;
  element3=[make_elem([n1 n2 n2+1],numu-1,1,1,numu);
	    make_elem([n1 n2+1 n1+1],numu-1,1,1,numu)];
	    
  element=[element1;element2;element3];
  numelem=rows(element);  

end

% plot the mesh
clf
plot_mesh(node,element,type)

% get the boundaries
b1=1:numu;
b2=numu:numu:numu*numv;
b3=[numu*numv:numu:numu*(numv+numv2-2),...
    numnode:-numu:offset+numu*numv];
b4=offset+numu*numv:-numu:offset+numu;
b5=offset+numu:-1:offset+1;
b6=[offset+1:numu:numnode-numu+1,numu*(numv+numv2-3)+1:-numu:1];

