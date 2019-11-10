knotVec     = [0 0 1 1];
controlPts  = [0 0; 1 0];
weights     = [1 1]; 

noCtrPts    = size(controlPts,1);
cp          = zeros(4,noCtrPts);
cp(1:2,:)   = controlPts';
cp(4,:)     = weights;

curve = nrbmak(cp,knotVec);

curve = nrbdegelev(curve,2); 

knots = [0.42 0.42 0.42 0.5 0.5 0.5 0.58 0.58 0.58];

curve     = nrbkntins(curve,knots);

uKnot     = curve.knots;


knotVec    = curve.knots;
p          = curve.order-1;
controlPts = curve.coefs(1:2,:)';
weights    = curve.coefs(4,:);
noCtrPts   = size(controlPts,1);
