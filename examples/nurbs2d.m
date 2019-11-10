addpath ../C_files/
addpath ../fem_util/

clear all

uKnot      = [0 0 0 0.5 0.5 1 1 1];
vKnot      = [0 0  1 1];

controlPts = [0 0 0; 0 0 1; 1 0 1;2 0 1;2 0 0;
              0 2 0; 0 2 1; 1 2 1;2 2 1;2 2 0];

p          = 2;
q          = 1;

noPtsX = 5;
noPtsY = 2;

weights    = ones(1,10)';

refineCount = 4;

hRefinement2d 

fileName = 'try1.eps';

%plotMesh3 (controlPts,weights, uKnot, vKnot,...
%                   p,q,20, 'r-', fileName);
               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


uKnot      = [0 0 0 1 2 3 4 4 4];
vKnot      = [0 0  1 1];

controlPts = [0 0 0;1 1 0;2 4 0;3 4 0;4 1 0;5 0 0;
              0 0 4;1 1 4;2 4 4;3 4 4;4 1 4;5 0 4];

p          = 2;
q          = 1;

noPtsX = 6;
noPtsY = 2;

weights    = ones(1,noPtsX*noPtsY)';

weights(3) = 0.5;
weights(4) = 0.5;
weights(9) = 0.5;
weights(10) = 0.5;

uKnot = uKnot/max(uKnot);

refineCount = 1;

hRefinement2d 

fileName = 'try1.eps';

plotMesh3 (controlPts,weights, uKnot, vKnot,...
                   p,q,40, 'r-', fileName);             
              





