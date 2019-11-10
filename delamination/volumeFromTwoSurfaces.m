function vol = volumeFromTwoSurfaces(surf1,surf2)

% make a volume from two surfaces

volumePts = zeros(4,surf.number(1),2,2);

volumePts(1:4,:,1,1) = surf1.coefs(:,:,1);
volumePts(1:4,:,2,1) = surf1.coefs(:,:,2);

volumePts(1:4,:,1,2) = surf2.coefs(:,:,1);
volumePts(1:4,:,2,2) = surf2.coefs(:,:,2);

uKnot = surface.knots{1};
vKnot = [0 0 1 1];
wKnot = [0 0 1 1];

vol   = nrbmak(volumePts,{uKnot vKnot wKnot});