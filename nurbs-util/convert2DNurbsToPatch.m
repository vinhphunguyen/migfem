function patch = convert2DNurbsToPatch(nurbs)

% convert a 2D NURBS object (created using nrbmak(controlPts,{uKnot
% vKnot});) to the usual data structure
% Vinh Phu Nguyen, March 2012
% nvinhphu@gmail.com

p          = nurbs.order(1)-1;
q          = nurbs.order(2)-1;
uKnot      = cell2mat(nurbs.knots(1));
vKnot      = cell2mat(nurbs.knots(2));
noPtsX     = length(uKnot)-p-1;
noPtsY     = length(vKnot)-q-1;
weights    = reshape(nurbs.coefs(4,:,:),noPtsX*noPtsY,1);

controlPts = [];

for iy=1:noPtsY
    controlPts = [controlPts; nurbs.coefs(1:2,:,iy)'];
end

controlPts(:,1)=controlPts(:,1)./weights;
controlPts(:,2)=controlPts(:,2)./weights;

patch = patch2D(uKnot,vKnot,p,q,controlPts,weights);