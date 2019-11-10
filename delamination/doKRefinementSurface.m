function solid = doKRefinementSurface(solid,pnew,qnew,refCount)

orders = solid.order - 1;
p      = orders(1);
q      = orders(2);

solid  = nrbdegelev(solid,[pnew-p qnew-q]); 

uKnot  = cell2mat(solid.knots(1));
vKnot  = cell2mat(solid.knots(2));

for i=1:refCount
    uKnotVectorU = unique(uKnot);
    uKnotVectorV = unique(vKnot);   
    newKnotsY = [];    
    % new knots along two directions (uniform)    
    newKnotsX = uKnotVectorU(1:end-1) + 0.5*diff(uKnotVectorU);
    %newKnotsY = uKnotVectorV(1:end-1) + 0.5*diff(uKnotVectorV);    
    newKnots  = {newKnotsX newKnotsY};    
    % h-refinement    
    solid     = nrbkntins(solid,newKnots);    
    uKnot     = cell2mat(solid.knots(1));
    vKnot     = cell2mat(solid.knots(2));
end