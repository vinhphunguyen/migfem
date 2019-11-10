% Data for the Scordelis-Lo roof problem
% Vinh Phu Nguyen, Feb 2013
% Delft University of Technology

addpath ../nurbs-geopdes/inst/

R      = 2540;
L      = 254;
t      = 12.7; % thickness 6.35
phi    = 0.1;  % radian

solid = nrbcylind(L,R,[],pi/2-phi,deg2rad(90));

%% refinement

% p-refinement to raise order

solid = nrbdegelev(solid,[1 2]);

uKnot      = cell2mat(solid.knots(1));
vKnot      = cell2mat(solid.knots(2));

refineCount = 1;

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

%% convert NURBS data back to our data structure for analysis

convert2DNurbsShell

%% Material properties

E  = 3102.75;
nu = 0.3;


%% Boundary condition




