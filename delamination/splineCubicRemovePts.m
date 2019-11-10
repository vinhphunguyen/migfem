addpath ../C_files/
addpath ../fem_util/
addpath ../examples/
addpath ../meshing/
addpath ../nurbs-util/

%%%%
% CUbic Bezier with offset distance larger than curvature radius!!!

clear all

% backtracking parameters
alpha = 0.1;
beta  = 0.7;

% initial curve

knotVec     = [0 0 0 1 1 1];
controlPts  = [0 0; -0.5 -4;3 -1];
weights     = [1 1 1]; % b-spline curves
p           = 2;

knotVecNew  = [0 0 0 0.3 0.6 0.8 1 1 1];
weightsNew  = [1 1 1 1 1 1]; 
pNew        = 2;

noCtrPts    = size(controlPts,1);
cp          = zeros(4,noCtrPts);
cp(1:2,:)   = controlPts';
cp(4,:)     = weights;

originalCurve = nrbmak(cp,knotVec);


figure
hold on
nrbctrlplot(originalCurve);
axis equal
axis off

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,'splinecurve.eps',opts)

% generating points on offset curve
thickness = -1.;

% computing offset points

noPts=100;
xi=linspace(0,1,noPts);
offsetPts = zeros(noPts,2);

for i=1:length(xi)
    xii = xi(i);
    [N dNdxi] = NURBS1DBasisDers(xii,p,knotVec,weights);    
    s   = findspan(noCtrPts-1,p,xii,knotVec);    
    pts = controlPts(s-p+1:s-p+1+p,:);
    x   = N    * pts;
    dx  = dNdxi* pts;
    dx  = dx/norm(dx);
    offsetPts(i,:) = x - [-dx(2) dx(1)]*thickness;
end

% remove intersecting offset points
offsetPts(22:72,:) = [];
noPts = size(offsetPts,1);
xi=linspace(0,1,noPts);

figure
hold on
nrbctrlplot(originalCurve);
axis equal
axis off
plot(offsetPts(:,1),offsetPts(:,2),'r*-','LineWidth',1.3);

%% gradient descent method 
% initial guess curve = a line

cpoints        = zeros(length(weightsNew),2);

cpoints(1,  :) = offsetPts(1,:);
cpoints(end,:) = offsetPts(end,:);

vec = cpoints(end,:) - cpoints(1,:);

dlam = 1/(length(weightsNew)-1);
for i=1:length(weightsNew)-2    
    cpoints(1+i,:) = cpoints(1,:)+ dlam*i*vec;    
end

cpoints0 = cpoints;

plot(cpoints0(:,1),cpoints0(:,2),'r-o',...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','g',...
    'MarkerSize',10);


samPts = computePoints(knotVecNew,weightsNew,cpoints,xi);

figure
hold on
axis equal
axis off
plot(offsetPts(:,1),offsetPts(:,2),'r*-','LineWidth',1.3);
plot(samPts(:,1),samPts(:,2),'bs-','LineWidth',1.2);

energy = getEnergy(offsetPts,samPts);

% disturbe the guess curve to compute numerically the gradient of the
% energy

h      = 1e-8;
eps    = 0.01;
eps1   = 0.01;
error  = 10;
error1 = 10;
dgamma = 0.0001;
gamma  = 0;
iiter  = 1;

while error > eps
    samPts   = computePoints(knotVecNew,weightsNew,cpoints,xi);
    f        = getEnergy(offsetPts,samPts);
    
    % determine the gradients numerically
    
    for i = 1:length(weightsNew)-2
        cpoints(1+i,:) = cpoints(1+i,:) + [h 0];
        samPts   = computePoints(knotVecNew,weightsNew,cpoints,xi);
        energy1x = getEnergy(offsetPts,samPts);
        cpoints(1+i,:) = cpoints(1+i,:) - [h 0] + [0 h];
        samPts   = computePoints(knotVecNew,weightsNew,cpoints,xi);
        energy1y = getEnergy(offsetPts,samPts);
        
        dir = 1/h*[energy1x-energy energy1y-energy];
        grad(i,:) = dir;
        cpoints(i+1,:) = cpoints0(i+1,:);
    end
    
    s = 1;
    for k=1:10
        cpnew = cpoints0;
        cpnew(2:end-1,:) = cpnew(2:end-1,:) - grad*s;
        samPts   = computePoints(knotVecNew,weightsNew,cpnew,xi);
        fnew     = getEnergy(offsetPts,samPts);
        if (fnew < f + s*alpha*(-grad)'*grad)
            break;
        else
            s = s*beta;
        end
    end
    
    cpoints0(2:end-1,:) = cpoints0(2:end-1,:) - grad*s;
    cpoints = cpoints0;
    
    samPts   = computePoints(knotVecNew,weightsNew,cpoints,xi);
    energy   = getEnergy(offsetPts,samPts);
    error    = energy;
    
    disp (sprintf(' %s %i %s %5.4e ', 'Iter',iiter, ':', error) );
    iiter = iiter + 1;
end

%% plot result

figure
hold on
plot(offsetPts(:,1),offsetPts(:,2),'r--','LineWidth',1.3);
%plot(samPts(:,1),   samPts(:,2),'cy--','LineWidth',1.3);

nrbctrlplot(originalCurve);
axis equal

cp        = zeros(4,size(cpoints,1));
cp(1:2,:) = cpoints0';
cp(4,:)   = weightsNew;
offsetCurve = nrbmak(cp,knotVecNew);
nrbctrlplot(offsetCurve);

%% 
% extruded surface

srf1 = nrbextrude(originalCurve, [0,0,1]);
srf2 = nrbextrude(offsetCurve, [0,0,1]);

figure
hold on
nrbctrlplot(srf1);
nrbctrlplot(srf2);

% make a volume from two surfaces

volumePts = zeros(4,noCtrPts,2,2);

volumePts(1:4,:,1,1) = srf1.coefs(:,:,1);
volumePts(1:4,:,2,1) = srf1.coefs(:,:,2);

volumePts(1:4,:,1,2) = srf2.coefs(:,:,1);
volumePts(1:4,:,2,2) = srf2.coefs(:,:,2);

uKnot = knotVec;
vKnot = [0 0 1 1];
wKnot = [0 0 1 1];

vol   = nrbmak(volumePts,{uKnot vKnot wKnot});

nrbplot(vol,[20 2 2])
axis off

%% continuum elements in between the curves

solidPts = zeros(4,noCtrPts,2);

solidPts(1:2,:,1) = controlPts';
solidPts(1:2,:,2) = cpoints';
solidPts(4,:)     = 1;

uKnot = knotVec;
vKnot = [0 0 1 1];

solid = nrbmak(solidPts,{uKnot vKnot});

figure 
hold on
nrbctrlplot(solid)
axis off

refineLevel = 3;
for i=1:refineLevel
    uKnotVectorU = unique(uKnot);
    uKnotVectorV = unique(vKnot);
    
    % new knots along two directions
    
    newKnotsX = uKnotVectorU(1:end-1) + 0.5*diff(uKnotVectorU);
    newKnotsY = uKnotVectorV(1:end-1) + 0.5*diff(uKnotVectorV);
    
    newKnots  = {newKnotsX []};
    solid     = nrbkntins(solid,newKnots);
    uKnot     = cell2mat(solid.knots(1));
    vKnot     = cell2mat(solid.knots(2));
end
%% 

convert2DNurbs

plotMesh (controlPts,weights,uKnot,vKnot,p,q,100,'r--','try.eps');







