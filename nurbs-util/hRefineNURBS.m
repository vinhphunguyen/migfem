function nurbs = hRefineNURBS(nurbs,refineCount)

% h-refinement for NURBS surfaces
% Vinh Phu Nguyen, March 2012
% nvinhphu@gmail.com

uKnot      = cell2mat(nurbs.knots(1));
vKnot      = cell2mat(nurbs.knots(2));

for i=1:refineCount
    uKnotVectorU = unique(uKnot);
    uKnotVectorV = unique(vKnot);
    
    % new knots along two directions (uniform)
    
    newKnotsX = uKnotVectorU(1:end-1) + 0.5*diff(uKnotVectorU);
    newKnotsY = uKnotVectorV(1:end-1) + 0.5*diff(uKnotVectorV);
    
    newKnots  = {newKnotsX newKnotsY};
    
    % h-refinement
    
    nurbs     = nrbkntins(nurbs,newKnots);    
    uKnot     = cell2mat(nurbs.knots(1));
    vKnot     = cell2mat(nurbs.knots(2));
end