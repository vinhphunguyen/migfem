function samPts = computePoints(knots,weights,cpoints,xi)

cp = zeros(4,size(cpoints,1));
cp(1:2,:) = cpoints';
cp(4,:)   = weights;
cp(1,:) = cp(1,:) .* weights;
cp(2,:) = cp(2,:) .* weights;

curve  = nrbmak(cp,knots);
[samPts] = nrbeval(curve,xi);
samPts = samPts(1:2,:)';