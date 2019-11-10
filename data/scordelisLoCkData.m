% Data for the Scordelis-Lo roof problem
% Quadratic NURBS for thickness direction
% Other directions: k-refinement using the NURBS toolbox
% Vinh Phu Nguyen, March 2012
% Delft University of Technology

addpath ../nurbs-geopdes/inst/

R      = 25;
L      = 25;
t      = 0.25;
phi    = 40;
Ri     = R-t;
Rii    = 0.5*(R+Ri);

%% control points

controlPts = zeros(4,3,3,3);

% first zKnot 

x1 = sqrt(R*R     + R*R*tan(deg2rad(20))^2);
x2 = sqrt(Rii*Rii + Rii*Rii*tan(deg2rad(20))^2);
x3 = sqrt(Ri*Ri   + Ri*Ri*tan(deg2rad(20))^2);

controlPts(1:3,1,1,1) = [R*cos(deg2rad(90-phi));R*sin(deg2rad(90-phi));0];
controlPts(1:3,2,1,1) = [x1*cos(deg2rad(70));x1*sin(deg2rad(70));0];
controlPts(1:3,3,1,1) = [0;R;0];

controlPts(1:3,1,2,1) = [Rii*cos(deg2rad(90-phi));Rii*sin(deg2rad(90-phi));0];
controlPts(1:3,2,2,1) = [x2*cos(deg2rad(70));x2*sin(deg2rad(70));0];
controlPts(1:3,3,2,1) = [0;Rii;0];

controlPts(1:3,1,3,1) = [Ri*cos(deg2rad(90-phi));Ri*sin(deg2rad(90-phi));0];
controlPts(1:3,2,3,1) = [x3*cos(deg2rad(70));x3*sin(deg2rad(70));0];
controlPts(1:3,3,3,1) = [0;Ri;0];

% second zKnot 
z = L/2;

controlPts(1:3,1,1,2) = [R*cos(deg2rad(90-phi));R*sin(deg2rad(90-phi));z];
controlPts(1:3,2,1,2) = [x1*cos(deg2rad(70));x1*sin(deg2rad(70));z];
controlPts(1:3,3,1,2) = [0;R;z];

controlPts(1:3,1,2,2) = [Rii*cos(deg2rad(90-phi));Rii*sin(deg2rad(90-phi));z];
controlPts(1:3,2,2,2) = [x2*cos(deg2rad(70));x2*sin(deg2rad(70));z];
controlPts(1:3,3,2,2) = [0;Rii;z];

controlPts(1:3,1,3,2) = [Ri*cos(deg2rad(90-phi));Ri*sin(deg2rad(90-phi));z];
controlPts(1:3,2,3,2) = [x3*cos(deg2rad(70));x3*sin(deg2rad(70));z];
controlPts(1:3,3,3,2) = [0;Ri;z];

% third zKnot
z = L;

controlPts(1:3,1,1,3) = [R*cos(deg2rad(90-phi));R*sin(deg2rad(90-phi));z];
controlPts(1:3,2,1,3) = [x1*cos(deg2rad(70));x1*sin(deg2rad(70));z];
controlPts(1:3,3,1,3) = [0;R;z];

controlPts(1:3,1,2,3) = [Rii*cos(deg2rad(90-phi));Rii*sin(deg2rad(90-phi));z];
controlPts(1:3,2,2,3) = [x2*cos(deg2rad(70));x2*sin(deg2rad(70));z];
controlPts(1:3,3,2,3) = [0;Rii;z];

controlPts(1:3,1,3,3) = [Ri*cos(deg2rad(90-phi));Ri*sin(deg2rad(90-phi));z];
controlPts(1:3,2,3,3) = [x3*cos(deg2rad(70));x3*sin(deg2rad(70));z];
controlPts(1:3,3,3,3) = [0;Ri;z];

controlPts(4,:,:,:)   = 1;

fac                   = cos(deg2rad(phi/2));

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
solid     = nrbkntins(solid,newKnots);
    
%% refinement

% p-refinement to raise order

%solid = nrbdegelev(solid,[3 0 3]); 

refineCount = 4;

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
nrbplot(solid,[80 1 80])
view([0 90])

crv = nrbcirc(R,[],deg2rad(0),deg2rad(180));
nrbplot(crv,80);
crv = nrbcirc(R-t,[],deg2rad(0),deg2rad(180));
nrbplot(crv,80);

%%%%%%%%
%% convert NURBS data back to our data structure for analysis

convert3DNurbs

%%%

buildVisualization3dMesh
figure
hold on
plot_mesh(node,elementV,'B8','b.-',1.1);
view(3)
plot3(controlPts(:,1),controlPts(:,2),controlPts(:,3),'r*');

numNode = size(node,1);

% normal stresses
sigmaXX = zeros(numNode,1);
sigmaYY = zeros(numNode,1);
sigmaZZ = zeros(numNode,1);

% shear stresses
sigmaXY = zeros(numNode,1);
sigmaYZ = zeros(numNode,1);
sigmaZX = zeros(numNode,1);

% displacements
dispX   = zeros(numNode,1);
dispY   = zeros(numNode,1);
dispZ   = zeros(numNode,1);


VTKPostProcess3d(node,elementV,'B8','pinchedCylinderTemp',...
    [sigmaXX sigmaYY sigmaZZ sigmaXY sigmaYZ sigmaZX],[dispX dispY dispZ]);
