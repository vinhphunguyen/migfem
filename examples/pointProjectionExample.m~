addpath C_files/
addpath('~/code/xfem-efg-matlab/fem_util');

clear all

knots      = [0 0 0 1 1 1];
controlPts = [0 0; 0.3 1; 1 0];
p          = 2;
weights    = [1 1 1]; % NURBS curves


x  = NURBSinterpolation(0.33, p, knots, controlPts(:,1), weights)

u0 = 0.4;
%xi = pointProjection(u0,0.5,p,knots,[controlPts weights']);




