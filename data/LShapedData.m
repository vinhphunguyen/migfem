% NURBS data for L shaped sample
% NURBS data: control points, knots, basis orders and weights
% Vinh Phu Nguyen, Johns Hopkins University

p = 1;
q = 1;

uKnot = [0 0 0.5 1 1];
vKnot  = [0 0 1 1];

controlPts = [ -1 1; -1 -1; 1 -1;
               0 1;0 0; 1 0];


noPtsX =  3;
noPtsY =  2;

weights = ones(1,noPtsX*noPtsY)';
