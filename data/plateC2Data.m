% rectangular plate in tension

controlPts          = zeros(4,2,2);

controlPts(1:2,1,1) = [0;0];
controlPts(1:2,2,1) = [1;0];
controlPts(1:2,1,2) = [0;1];
controlPts(1:2,2,2) = [1;1];

controlPts(4,:,:)   = 1;

uKnot = [0 0 1 1];
vKnot = [0 0 1 1];

solid = nrbmak(controlPts,{uKnot vKnot});

% evaluate order

solid = nrbdegelev(solid,[2 2]);

% h-refinement (knot insertion)

uKnot     = cell2mat(solid.knots(1));
vKnot     = cell2mat(solid.knots(2));

refineCountX = 3;
for i=1:refineCountX
    uKnotVectorU = unique(uKnot);
    uKnotVectorV = unique(vKnot);
    newKnotsX = uKnotVectorU(1:end-1) + 0.5*diff(uKnotVectorU);
    newKnotsY = uKnotVectorV(1:end-1) + 0.5*diff(uKnotVectorV);
    newKnots  = {newKnotsX newKnotsY};
    solid     = nrbkntins(solid,newKnots);
    uKnot      = cell2mat(solid.knots(1));
    vKnot      = cell2mat(solid.knots(2));
end


