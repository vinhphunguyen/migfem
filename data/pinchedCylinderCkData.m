% Data for the pinched cylinder problem
% Quadratic NURBS for thickness direction
% Other directions: k-refinement using the NURBS toolbox
% Vinh Phu Nguyen, March 2012
% Delft University of Technology

addpath ../nurbs-geopdes/inst/
addpath ../post-processing/
addpath ../

R0     = 300;
L      = 300;
t      = 3;
R      = R0 +t/2;
Ri     = R0 -t/2;
Rii    = 0.5*(R+Ri);

%% control points

controlPts = zeros(4,3,3,3);

% first zKnot 

controlPts(1:3,1,1,1) = [R;0;0];
controlPts(1:3,2,1,1) = [R;R;0];
controlPts(1:3,3,1,1) = [0;R;0];

controlPts(1:3,1,2,1) = [Rii;0;0];
controlPts(1:3,2,2,1) = [Rii; Rii;0];
controlPts(1:3,3,2,1) = [0;Rii;0];

controlPts(1:3,1,3,1) = [Ri;0;0];
controlPts(1:3,2,3,1) = [Ri; Ri;0];
controlPts(1:3,3,3,1) = [0;Ri;0];

% second zKnot 
z = L/2;
controlPts(1:3,1,1,2) = [R;0;z];
controlPts(1:3,2,1,2) = [R; R;z];
controlPts(1:3,3,1,2) = [0;R;z];

controlPts(1:3,1,2,2) = [Rii;0;z];
controlPts(1:3,2,2,2) = [Rii;Rii;z];
controlPts(1:3,3,2,2) = [0;Rii;z];

controlPts(1:3,1,3,2) = [Ri;0;z];
controlPts(1:3,2,3,2) = [Ri;Ri;z];
controlPts(1:3,3,3,2) = [0;Ri;z];

% third zKnot
z = L;
controlPts(1:3,1,1,3) = [R;0;z];
controlPts(1:3,2,1,3) = [R; R;z];
controlPts(1:3,3,1,3) = [0;R;z];

controlPts(1:3,1,2,3) = [Rii;0;z];
controlPts(1:3,2,2,3) = [Rii;Rii;z];
controlPts(1:3,3,2,3) = [0;Rii;z];

controlPts(1:3,1,3,3) = [Ri;0;z];
controlPts(1:3,2,3,3) = [Ri; Ri;z];
controlPts(1:3,3,3,3) = [0;Ri;z];

controlPts(4,:,:,:)   = 1;

fac                   = 1/sqrt(2);

controlPts(4,2,1,1) = fac;
controlPts(4,2,2,1) = fac;
controlPts(4,2,3,1) = fac;

controlPts(4,2,1,2) = fac;
controlPts(4,2,2,2) = fac;
controlPts(4,2,3,2) = fac;

controlPts(4,2,1,3) = fac;
controlPts(4,2,2,3) = fac;
controlPts(4,2,3,3) = fac;

% homogenous coordinates (x*w,y*w,z*w)

controlPts(1:3,2,1,1) = controlPts(1:3,2,1,1)*fac;
controlPts(1:3,2,2,1) = controlPts(1:3,2,2,1)*fac;
controlPts(1:3,2,3,1) = controlPts(1:3,2,3,1)*fac;
controlPts(1:3,2,1,2) = controlPts(1:3,2,1,2)*fac;
controlPts(1:3,2,2,2) = controlPts(1:3,2,2,2)*fac;
controlPts(1:3,2,3,2) = controlPts(1:3,2,3,2)*fac;
controlPts(1:3,2,1,3) = controlPts(1:3,2,1,3)*fac;
controlPts(1:3,2,2,3) = controlPts(1:3,2,2,3)*fac;
controlPts(1:3,2,3,3) = controlPts(1:3,2,3,3)*fac;


%% knot vectors 1x1x1 mesh

uKnot = [0 0 0 1 1 1];
vKnot = [0 0 0 1 1 1];
wKnot = [0 0 0 1 1 1];

%% build NURBS object

solid = nrbmak(controlPts,{uKnot vKnot wKnot});

newKnots  = {[] [0.5] []};
%solid     = nrbkntins(solid,newKnots);

%% k-refinement

% p-refinement to raise order

%solid = nrbdegelev(solid,[1 0 1]); 

% then h-refinement 

refineCount = 0;

for i=1:refineCount
    uKnotVectorU = unique(uKnot);
    uKnotVectorW = unique(wKnot);
    
    % new knots along two directions (uniform)
    
    newKnotsX = uKnotVectorU(1:end-1) + 0.5*diff(uKnotVectorU);
    newKnotsZ = uKnotVectorW(1:end-1) + 0.5*diff(uKnotVectorW);
    newKnots  = {newKnotsX [] newKnotsZ};
    
    % h-refinement
    
    solid     = nrbkntins(solid,newKnots);
    
    uKnot      = cell2mat(solid.knots(1));
    wKnot      = cell2mat(solid.knots(3));
end

figure 
hold on
%nrbplot(solid,[40 1 40])
%view([0 90])

nrbkntplot (solid)
view(3)

% crv = nrbcirc(R,[],deg2rad(0),deg2rad(180));
% nrbplot(crv,80);
% crv = nrbcirc(R-t,[],deg2rad(0),deg2rad(180));
% nrbplot(crv,80);

%%%%%%%%
%% convert NURBS data back to our data structure for analysis

convert3DNurbs

plotMesh3 (controlPts,weights, uKnot,vKnot,wKnot,...
                   p,q,r,100,'r-','try.eps')

%%%

buildVisualization3dMesh
figure
hold on
plot_mesh(node,elementV,'B8','b.-',1.1);
view(3)
view([0 90])
% plot3(controlPts(:,1),controlPts(:,2),controlPts(:,3),'r*');
% 
% numNode = size(node,1);
% 
% % normal stresses
% sigmaXX = zeros(numNode,1);
% sigmaYY = zeros(numNode,1);
% sigmaZZ = zeros(numNode,1);
% 
% % shear stresses
% sigmaXY = zeros(numNode,1);
% sigmaYZ = zeros(numNode,1);
% sigmaZX = zeros(numNode,1);
% 
% % displacements
% dispX   = zeros(numNode,1);
% dispY   = zeros(numNode,1);
% dispZ   = zeros(numNode,1);
% 
% 
% VTKPostProcess3d(node,elementV,'B8','pinchedCylinderTemp',...
%     [sigmaXX sigmaYY sigmaZZ sigmaXY sigmaYZ sigmaZX],[dispX dispY dispZ]);
