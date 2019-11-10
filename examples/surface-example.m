
controlPts = zeros(4,4,4);

controlPts(1:3,1,1) = [0;0;0];
controlPts(1:3,2,1) = [1;0;1.1];
controlPts(1:3,3,1) = [2;0;1.1];
controlPts(1:3,4,1) = [3;0;0];

controlPts(1:3,1,2) = [0;1;2];
controlPts(1:3,2,2) = [1;1;2];
controlPts(1:3,3,2) = [2;1;2];
controlPts(1:3,4,2) = [3;1;2];

controlPts(1:3,1,3) = [0;2;2];
controlPts(1:3,2,3) = [1;2;2];
controlPts(1:3,3,3) = [2;2;2];
controlPts(1:3,4,3) = [3;2;2];

controlPts(1:3,1,4) = [0;3;0];
controlPts(1:3,2,4) = [1;3;1.1];
controlPts(1:3,3,4) = [2;3;1.1];
controlPts(1:3,4,4) = [3;3;0];

controlPts(4,:) = 1;

uKnot = [0 0 0 0.5 1 1 1];
vKnot = [0 0 0 0.5 1 1 1];

surf = nrbmak(controlPts,{uKnot vKnot});


figure
hold on
nrbctrlplot(surf);
axis equal
axis off

