% convert a 2D NURBS object (created using nrbmak(controlPts,{uKnot
% vKnot});) to the usual data structure
% Vinh Phu Nguyen, March 2012
% nvinhphu@gmail.com

p          = solid.order(1)-1;
q          = solid.order(2)-1;
uKnot      = cell2mat(solid.knots(1));
vKnot      = cell2mat(solid.knots(2));
noPtsX     = length(uKnot)-p-1;
noPtsY     = length(vKnot)-q-1;
weights    = reshape(solid.coefs(4,:,:),noPtsX*noPtsY,1);

controlPts = [];

for iy=1:noPtsY
    controlPts = [controlPts; solid.coefs(1:3,:,iy)'];
end

% our controlPts only stores (x,y,z) not (w*x,w*y,w*z)

controlPts(:,1) = controlPts(:,1)./weights;
controlPts(:,2) = controlPts(:,2)./weights;
controlPts(:,3) = controlPts(:,3)./weights; % for shell problems, z-coord exists
