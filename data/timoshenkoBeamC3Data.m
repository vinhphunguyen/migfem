% Timoshenko beam 48x12 
% Discretised with C3 (p=q=4) NURBS elements
% Control points obtained automatically from a structured
% Q4 mesh. 
% Vinh Phu Nguyen
% Delft University of Technology, The Netherlands

L = 48;
D = 12;

% no of control points along X and Y directions

noPtsX = 6;
noPtsY = 6;

gcoord     = meshRectangularCoord(L,D,noPtsX-1,noPtsY-1);
gcoord(:,2)=gcoord(:,2)-6;
controlPts = gcoord;

% basis orders

p = 4;
q = 4;

% weights

weights = ones(1,noPtsX*noPtsY)';

% knot vectors

knotUTemp = linspace(0,1,noPtsX-p+1);
knotVTemp = linspace(0,1,noPtsY-q+1);

uKnot = [0 0 0 0 knotUTemp 1 1 1 1];
vKnot = [0 0 0 0 knotVTemp 1 1 1 1];

noCtrPts   = noPtsX * noPtsY;
noDofs     = noCtrPts * 2;

% sometimes h-refinement process gives NAN new control pts
% simply remove them with the following 

controlPts  = controlPts(1:noCtrPts,:);
weights     = weights(1:noCtrPts);


