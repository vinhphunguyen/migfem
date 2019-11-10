function [R,dRdx,J] = getShapeGrads2D(pt,xiE,etaE,pts)

global p q uKnot vKnot weights

% compute coords in parameter space
Xi      = parent2ParametricSpace(xiE,pt(1));
Eta     = parent2ParametricSpace(etaE,pt(2));
J2      = jacobianPaPaMapping(xiE,etaE);

% compute derivatives of shape functions
[R dRdxi dRdeta] = NURBS2DBasisDers([Xi; Eta],p,q,uKnot,vKnot,weights');

% compute the jacobian of physical and parameter domain mapping
% then the derivative w.r.t spatial physical coordinates

jacob      = [dRdxi; dRdeta] * pts;
J1         = det(jacob);
invJacob   = inv(jacob);
dRdx       = invJacob * [dRdxi; dRdeta];
J          = J1 * J2;