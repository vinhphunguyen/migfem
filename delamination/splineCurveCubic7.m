addpath ../C_files/
addpath ../fem_util/
addpath ../examples/
addpath ../meshing/
addpath ../nurbs-geopdes/inst/

clear all

% backtracking parameters
alpha = 0.1;
beta  = 0.7;

% initial curve

knotVec     = [0 0 0 0 1 2 3 4 4 4 4];
knotVec     = knotVec/max(knotVec);
controlPts  = [0.5 0;0.5 1;2 0.5;3 -1;2 -1;0.5 -0.9; 0.5 0];
weights     = [1 1 1 1 1 1 1]; % b-spline curves
p           = 3;

noCtrPts    = size(controlPts,1);
cp          = zeros(4,noCtrPts);
cp(1:2,:)   = controlPts';
cp(4,:)     = weights;

originalCurve = nrbmak(cp,knotVec);

figure
hold on
nrbctrlplot(originalCurve);
axis equal


opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,'splinecurve.eps',opts)

% generating points on offset curve
thickness = 0.35;

% computing offset points
%refinement = 0;
%generateIGA1DMesh;
%noElems   = size(elConn,1);    % no of elements
%noGPs     = 30;
%noPts     = noElems*noGPs;
%[W,Q]     = quadrature(  noGPs, 'GAUSS', 1 ); 


noPts=120;
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
    offsetPts(i,:) = x + [-dx(2) dx(1)]*thickness;
end

%plot(offsetPts(:,1),offsetPts(:,2),'r*-','LineWidth',1.3);

%% gradient descent method 
% initial guess curve = a line

cpoints        = controlPts(:,1:2);

cpoints(1,  :) = offsetPts(1,:);
cpoints(end,:) = offsetPts(end,:);

vec = cpoints(end,:) - cpoints(1,:);

dlam = 1/(noCtrPts-1);
for i=1:noCtrPts-2    
    cpoints(1+i,:) = cpoints(1,:)+ dlam*i*vec;    
end

cpoints0 = cpoints;

plot(cpoints0(:,1),cpoints0(:,2),'r-o',...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','g',...
    'MarkerSize',10);


samPts = computePoints(knotVec,weights,cpoints,xi);

plot(samPts(:,1),samPts(:,2),'bs-','LineWidth',1.8);

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
    samPts   = computePoints(knotVec,weights,cpoints,xi);
    f        = getEnergy(offsetPts,samPts);
    
    % determine the gradients numerically
    
    for i = 1:noCtrPts-2
        cpoints(1+i,:) = cpoints(1+i,:) + [h 0];
        samPts   = computePoints(knotVec,weights,cpoints,xi);
        energy1x = getEnergy(offsetPts,samPts);
        cpoints(1+i,:) = cpoints(1+i,:) - [h 0] + [0 h];
        samPts   = computePoints(knotVec,weights,cpoints,xi);
        energy1y = getEnergy(offsetPts,samPts);
        
        dir = 1/h*[energy1x-energy energy1y-energy];
        grad(i,:) = dir;
        cpoints(i+1,:) = cpoints0(i+1,:);
    end
    
    s = 1;
    for k=1:10
        cpnew = cpoints0;
        cpnew(2:end-1,:) = cpnew(2:end-1,:) - grad*s;
        samPts   = computePoints(knotVec,weights,cpnew,xi);
        fnew     = getEnergy(offsetPts,samPts);
        if (fnew < f + s*alpha*(-grad)'*grad)
            break;
        else
            s = s*beta;
        end
    end
    
    cpoints0(2:end-1,:) = cpoints0(2:end-1,:) - grad*s;
    cpoints = cpoints0;
    
    samPts   = computePoints(knotVec,weights,cpoints,xi);
    energy   = getEnergy(offsetPts,samPts);
    error    = energy;
    
    disp (sprintf(' %s %i %s %5.4e ', 'Iter',iiter, ':', error) );
    iiter = iiter + 1;
end

%% plot result

figure
hold on
%plot(offsetPts(:,1),offsetPts(:,2),'r--','LineWidth',1.3);
%plot(samPts(:,1),   samPts(:,2),'cys--','LineWidth',1.3);

nrbplot(originalCurve,120);
axis equal

cp        = zeros(4,size(controlPts,1));
cp(1:2,:) = cpoints0';
cp(4,:)   = weights;
offsetCurve = nrbmak(cp,knotVec);
nrbplot(offsetCurve,120);

%% 
% extruded surface

solid = nrbextrude(originalCurve, [0,1]);
srf2 = nrbextrude(offsetCurve, [0,1]);

figure
hold on
nrbctrlplot(srf1);
nrbctrlplot(srf2);

convert2DNurbs
plotMesh (controlPts,weights,uKnot,vKnot,p,q,100,'b-','try.eps');


