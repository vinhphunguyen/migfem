function {node,element}=mirrormesh(node,element,p1,p2)

% {NODE,ELEMENT}=MIRRORMESH(NODE,ELEMENT,PT1,PT2)
%    
%

% mirror nodes and keep a look-up table
N=rows(node);
luTable=0*node;
TOL=.0001;

mirrorDir=[p1(2)-p2(2),p2(1)-p1(1)];

for n=1:N
  
  % find distance to mirror line  
  A=[
