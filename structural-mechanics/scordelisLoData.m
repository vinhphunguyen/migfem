% Data for the Scordelis-Lo roof problem
% Vinh Phu Nguyen, Feb 2013
% Delft University of Technology

addpath ../nurbs-geopdes/inst/

R      = 25;
L      = 25;
phi    = 40;

%% control points

controlPts = zeros(4,3,3);

% first zKnot 

x1 = sqrt(R*R + R*R*tan(deg2rad(20))^2);

controlPts(1:3,1,1) = [R*cos(deg2rad(90-phi));R*sin(deg2rad(90-phi));0];
controlPts(1:3,2,1) = [x1*cos(deg2rad(70));x1*sin(deg2rad(70));0];
controlPts(1:3,3,1) = [0;R;0];

% second zKnot 
z = L/2;

controlPts(1:3,1,2) = [R*cos(deg2rad(90-phi));R*sin(deg2rad(90-phi));z];
controlPts(1:3,2,2) = [x1*cos(deg2rad(70));x1*sin(deg2rad(70));z];
controlPts(1:3,3,2) = [0;R;z];

% third zKnot
z = L;

controlPts(1:3,1,3) = [R*cos(deg2rad(90-phi));R*sin(deg2rad(90-phi));z];
controlPts(1:3,2,3) = [x1*cos(deg2rad(70));x1*sin(deg2rad(70));z];
controlPts(1:3,3,3) = [0;R;z];

controlPts(4,:,:)   = 1; % weights

fac                 = cos(deg2rad(phi/2));

controlPts(4,2,1) = fac;
controlPts(4,2,2) = fac;
controlPts(4,2,3) = fac;

% homogenous coordinates (x*w,y*w,z*w)

controlPts(1:3,2,1) = controlPts(1:3,2,1)*fac;
controlPts(1:3,2,2) = controlPts(1:3,2,2)*fac;
controlPts(1:3,2,3) = controlPts(1:3,2,3)*fac;


%% knot vectors 1x1 mesh

uKnot = [0 0 0 1 1 1];
vKnot = [0 0 0 1 1 1];


%% build NURBS object

solid = nrbmak(controlPts,{uKnot vKnot});
    
%% refinement

% p-refinement to raise order

solid = nrbdegelev(solid,[1 1]); 

refineCount = 4;

for i=1:refineCount
    uKnotVectorU = unique(uKnot);
    uKnotVectorW = unique(vKnot);
    
    % new knots along two directions (uniform)
    
    newKnotsX = uKnotVectorU(1:end-1) + 0.5*diff(uKnotVectorU);
    newKnotsZ = uKnotVectorW(1:end-1) + 0.5*diff(uKnotVectorW);
    newKnots  = {newKnotsX newKnotsZ};
    
    % h-refinement
    
    solid     = nrbkntins(solid,newKnots);
    
    uKnot      = cell2mat(solid.knots(1));
    vKnot      = cell2mat(solid.knots(2));
end

figure 
hold on
nrbplot(solid,[40 40])
view(3)

crv = nrbcirc(R,[],deg2rad(0),deg2rad(180));
nrbplot(crv,80);
view([0 90])

%% convert NURBS data back to our data structure for analysis

convert2DNurbsShell

%% Material properties

E  = 4.32e8;
nu = 0.0;
t  = 0.25; % thickness

%% Boundary condition

g  = -90;
