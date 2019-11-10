function [mesh]=buildGrid3D(Lx,Ly,Lz,noX,noY,noZ)

controlPts          = zeros(4,2,2);

controlPts(1:2,1,1) = [0;0];
controlPts(1:2,2,1) = [Lx;0];
controlPts(1:2,1,2) = [0;Ly];
controlPts(1:2,2,2) = [Lx;Ly];

controlPts(4,:,:)   = 1;

uKnot = [0 0 1 1];
vKnot = [0 0 1 1];

nurbs  = nrbmak(controlPts,{uKnot vKnot});

% extrusion to have 3D objects
nurbs  = nrbextrude(nurbs, [0,0,Lz]);
nurbs  = doKRefinementSolid(nurbs,  [1 1 1], [noX noY noZ]);
igaMesh  = buildIGA3DMesh(nurbs);
mesh   = buildVisualizationMesh3D(nurbs);
mesh.noElemsX = igaMesh.noElemsU;
mesh.noElemsY = igaMesh.noElemsV;
mesh.noElemsZ = igaMesh.noElemsW;
