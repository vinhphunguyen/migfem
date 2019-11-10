
a = 0.3; % inner radius
b = 0.5; % outer radius

% quadratic NURBS

p = 2;
q = 2;

uKnot = [0 0 0 1 1 1];
vKnot = [0 0 0 1 1 1];

controlPts=[a 0;a a;0 a;
            0.5*(a+b) 0; 0.5*(a+b) 0.5*(a+b);0 0.5*(a+b); 
            b 0; b b; 0 b];

noPtsX = length(uKnot)-p-1;
noPtsY = length(vKnot)-q-1;

weights = ones(1,noPtsX*noPtsY)';

weights(2)=1/sqrt(2);
weights(5)=1/sqrt(2);
weights(8)=1/sqrt(2);

% refineCount=1;
% hRefinement2d 
% plotMesh (controlPts,weights,uKnot,vKnot,p,q,50,'r--','try.eps');




