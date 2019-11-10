% convert GeoPDEs 3D NURBS object to the data structures
% used in IGA code.
% VP Nguyen, 2012.

p          = solid.order(1)-1;
q          = solid.order(2)-1;
r          = solid.order(3)-1;
uKnot      = cell2mat(solid.knots(1));
vKnot      = cell2mat(solid.knots(2));
wKnot      = cell2mat(solid.knots(3));
noPtsX     = length(uKnot)-p-1;
noPtsY     = length(vKnot)-q-1;
noPtsZ     = length(wKnot)-r-1;
weights    = reshape(solid.coefs(4,:,:),noPtsX*noPtsY*noPtsZ,1);

controlPts = [];

for iz=1:noPtsZ
    for iy=1:noPtsY
        controlPts = [controlPts; solid.coefs(1:3,:,iy,iz)'];
    end
end

% our controlPts only stores (x,y,z) not (w*x,w*y,w*z)

controlPts(:,1) = controlPts(:,1)./weights;
controlPts(:,2) = controlPts(:,2)./weights;
controlPts(:,3) = controlPts(:,3)./weights;
