%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [controlPoints, weights] = proj2nurbs(projcoord)
%--------------------------------------------------------------
%function [nob, controlPoints, weightVector] = proj2nurbs(projcoord)
% transform projective coordinates into NURBS data
%INPUT:
% projcoord    : matrix with projective coordinates
%OUTPUT:
% nob          : # of basis function = # control points / weights
% controlPoints: vector of control points (1 per row)
% weightVector : column vector of weights
%--------------------------------------------------------------

dimension     = size(projcoord,2);
weights       = projcoord(:,dimension);
controlPoints = projcoord(:,1:dimension-1);

for i=1:size(weights,1)
    controlPoints(i,:) = controlPoints(i,:)* 1/(weights(i));
end
end