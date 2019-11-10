function solid = surfaceFromTwoCurves(curve1,curve2)

noCtrPts = curve1.number;

solidPts = zeros(4,noCtrPts,2);

solidPts(1:2,:,1) = curve1.coefs(1:2,:);
solidPts(1:2,:,2) = curve2.coefs(1:2,:);
solidPts(4,:)     = 1;

uKnot = curve1.knots;
vKnot = [0 0 1 1];

solid = nrbmak(solidPts,{uKnot vKnot});
