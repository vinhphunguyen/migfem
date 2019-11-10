% rectangular plate in tension


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

controlPts=[0 0; 0 0.5; 0 1;
            0.5 0;0.5 0.5;0.5 1;
            1 0; 1 0.5; 1 1];


% knot vectors

uKnot = [0 0 0 1 1 1];
vKnot = [0 0 0 1 1 1];

p = 2;
q = 2;

noPtsX = 3;
noPtsY = 3;

weights = ones(1,noPtsX*noPtsY)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

controlPts = [

         0         0;
         0.5    0.;
         1.0000 0;
    0.         0.5;    
    0.500    0.500;
    1    0.5;
    0        1;
    0.5   1;
    1.0000    1.0000];

uKnot = [0 0 0 1 1 1];
vKnot = [0 0 0 1 1 1];

p = 2;
q = 2;

noPtsX = 3;
noPtsY = 3;

weights = ones(1,noPtsX*noPtsY)';

