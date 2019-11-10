% NURBS data for L shaped sample
% NURBS data: control points, knots, basis orders and weights
% Vinh Phu Nguyen, Johns Hopkins University

p = 2;
q = 2;

uKnot = [0 0 0 0.5 1 1 1];
vKnot  = [0 0 0 1 1 1];
weights = [1 1 1 1 1 1 1 1 1 1 1 1]';

a = 50;

controlPts = [ 2*a 0;0 0; 0 0; 0 2*a;
               2*a 0.5*a; a 0.5*a;0.25*a a;0.25*a 2*a;
               2*a a; a a; a a; a 2*a];


noPtsX =  4;
noPtsY =  3;



