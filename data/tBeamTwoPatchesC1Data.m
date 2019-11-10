% Timoshenko beam 48x12 
% The beam is constructed by two patches.
% Vinh Phu Nguyen
% Delft University of Technology, The Netherlands

L = 48;
D = 12;
a = 24;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% patch1

% no of control points along X and Y directions

noPtsX = 20;
noPtsY = 8;

gcoord      = meshRectangularCoord(L-a,D,noPtsX-1,noPtsY-1);
gcoord(:,2) = gcoord(:,2)-6;
controlPts1 = gcoord;

% basis orders

p = 2;
q = 2;

% weights

weights = ones(1,noPtsX*noPtsY)';

% knot vectors

knotUTemp = linspace(0,1,noPtsX-p+1);
knotVTemp = linspace(0,1,noPtsY-q+1);

uKnot = [0 0 knotUTemp 1 1];
vKnot = [0 0 knotVTemp 1 1];

% build a patch object

patch1 = patch2D(uKnot,vKnot,p,q,controlPts1,weights);

% build mesh data

chan  = zeros(noPtsY,noPtsX);

count = 1;

for i=1:noPtsY
    for j=1:noPtsX
        chan(i,j) = count;
        count = count + 1;
    end
end

generateMesh(patch1,chan,chan);

% find control points on the interface between two patches

index = find(controlPts1(:,1)==L-a);

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% patch2

% no of control points along X and Y directions

noPtsX = 6;
noPtsY = 8;

gcoord      = meshRectangularCoord(a,D,noPtsX-1,noPtsY-1);
gcoord(:,1) = gcoord(:,1)+L-a;
gcoord(:,2) = gcoord(:,2)-6;
controlPts2 = gcoord;

% basis orders

p = 2;
q = 2;

% weights

weights = ones(1,noPtsX*noPtsY)';

% knot vectors

knotUTemp = linspace(0,1,noPtsX-p+1);
knotVTemp = linspace(0,1,noPtsY-q+1);

uKnot = [0 0 knotUTemp 1 1];
vKnot = [0 0 knotVTemp 1 1];

% build a patch object

patch2 = patch2D(uKnot,vKnot,p,q,controlPts2,weights);

% build mesh data

chan      = zeros(noPtsY,noPtsX);
chan(:,1) = index;
nodeId    = index(end); % maximum number of node of patch1

count = 1;

for i=1:noPtsY
    for j=2:noPtsX
        chan(i,j) = count + nodeId;
        count = count + 1;
    end
end

chanLocal  = zeros(noPtsY,noPtsX);

count = 1;

for i=1:noPtsY
    for j=1:noPtsX
        chanLocal(i,j) = count;
        count = count + 1;
    end
end

generateMesh(patch2,chan,chanLocal);

%%

noPatches  = 2;
patches(1) = patch1;
patches(2) = patch2;


