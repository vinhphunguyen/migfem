

a=10;
b=20;
r=3;

uKnot = [0 0 0 1 1 2 2 3 3 4 4 4];
vKnot = [0 0 1 1];

p     = 2;
q     = 1;

controlPts = [0 0;0 0.5*b-r;0 0.5*b-r; r 0.5*b-r;r 0.5*b;r 0.5*b+r;0 0.5*b+r;0 0.5*b+r;0 b;
              a 0; a b/8; a 2*b/8; a 3*b/8; a 4*b/8; a 5*b/8; a 6*b/8; a 7*b/8;a b];

noPtsX = length(uKnot)-p-1;
noPtsY = length(vKnot)-q-1;

uKnot = uKnot/max(uKnot);

weights = ones(1,noPtsX*noPtsY)';
weights([4,6]) = 1/sqrt(2);

% refineCount =1;
% hRefinement2d 

plotMesh (controlPts,weights,uKnot,vKnot,p,q,50,'r-','try.eps');