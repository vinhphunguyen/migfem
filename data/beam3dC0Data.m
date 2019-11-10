% NURBS data for a 3D beam 
% NURBS data: control points, knots, basis orders and weights
% Vinh Phu Nguyen, Johns Hopkins University


% control points
% in order not to hard wire the control points, use
% a structured tri-linear brick mesh

a = 10;
b = 4;
c = 2;

noPtsX = 21;
noPtsY = 7;
noPtsZ = 3;

[controlPts,elementVV]=makeB8mesh(a,b,c,noPtsX,noPtsY,noPtsZ);

% knot vectors

knotUTemp = linspace(0,1,noPtsX);
knotVTemp = linspace(0,1,noPtsY);
knotWTemp = linspace(0,1,noPtsZ);

uKnot = [0 knotUTemp 1];
vKnot = [0 knotVTemp 1];
wKnot = [0 knotWTemp 1];

% basis order

p = 1;
q = 1;
r = 1;

% weights

weights = ones(1,noPtsX*noPtsY*noPtsZ)';




          
