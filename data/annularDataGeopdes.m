a = 0.3; % inner radius
b = 0.5; % outer radius

res = 200; % resolution for plotting NURBS


uKnot = [0 0 1 1];
vKnot = [0 0 0 1 1 1];

controlPts          = zeros(4,2,3);

controlPts(1:2,1,1) = [a;0];
controlPts(1:2,2,1) = [b;0;];


controlPts(1:2,1,2) = [a;a];
controlPts(1:2,2,2) = [b;b];

controlPts(1:2,1,3) = [0;a];
controlPts(1:2,2,3) = [0;b];

controlPts(4,:,:)   = 1;

fac                 = 1/sqrt(2);

controlPts(4,1,2) = fac;
controlPts(4,2,2) = fac;

% homogenous coordinates (x*w,y*w,z*w)

controlPts(1:2,1,2) = controlPts(1:2,1,2)*fac;
controlPts(1:2,2,2) = controlPts(1:2,2,2)*fac;

%% build NURBS object

solid = nrbmak(controlPts,{uKnot vKnot});

figure 
hold on
nrbplot(solid,[40 40])

%% 

refineLevel = 3;
for i=1:refineLevel
    uKnotVectorU = unique(uKnot);
    uKnotVectorV = unique(vKnot);
    
    % new knots along two directions
    
    newKnotsX = uKnotVectorU(1:end-1) + 0.5*diff(uKnotVectorU);
    newKnotsY = uKnotVectorV(1:end-1) + 0.5*diff(uKnotVectorV);
    
    newKnots  = {newKnotsX newKnotsY};
    solid     = nrbkntins(solid,newKnots);
    uKnot      = cell2mat(solid.knots(1));
    vKnot      = cell2mat(solid.knots(2));
end

%% 

convert2DNurbs


plotMesh (controlPts,weights,uKnot,vKnot,p,q,res,'r--','try.eps');


noCtrPts       = noPtsX * noPtsY;
noDofs         = noCtrPts * 2;

% find boundary nodes for boundary conditions

fixedXNodes  =  find(controlPts(:,1)==0);
fixedYNodes  =  find(controlPts(:,2)==0);
forcedNodes  =  1:noPtsX:noCtrPts;

% build connectivity ...

generateIGA2DMesh

% build boundary mesh for force vector computation

rightEdgeMesh  = zeros(noElemsV,q+1);

for i=1:noElemsV
    rightEdgeMesh(i,:) = forcedNodes(i:i+q);
end

neumannMesh.element = rightEdgeMesh;
neumannMesh.p       = q;
neumannMesh.knot    = vKnot;
neumannMesh.range   = elRangeV;


