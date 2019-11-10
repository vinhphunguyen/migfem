%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function projcoord = nurb2proj(nob, controlPoints, weights)
%--------------------------------------------------------------
%function projcoord = nurb2proj(nob, controlPoints, weights)
% transform NURBS data into projective coordinates
%INPUT:
% nob          : # of basis function = # control points / weights
% controlPoints: vector of control points (1 per row)
% weights :    : column vector of weights
%OUTPUT:
% projcoord    : matrix with projective coordinates
%--------------------------------------------------------------
projcoord = controlPoints;
for i=1:nob
    projcoord(i,:) = projcoord(i,:)*weights(i);
end
projcoord = [projcoord, weights];
end