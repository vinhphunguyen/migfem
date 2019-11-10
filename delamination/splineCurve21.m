addpath ../C_files/
addpath ../fem_util/
addpath ../examples/

clear all

% backtracking parameters
alpha = 0.1;  
beta  = 0.7;

knots      = [0 0 0 0.5 1 1 1];

uniqueKnot = unique(knots);

controlPts  = [0 0; 1 0.5; 2 -0.4; 3 0];

cp = zeros(4,4);
cp(1:2,:) = controlPts';
cp(4,:)  = 1;

originalCurve = nrbmak(cp,knots);

p          = 2;

weights   = [1 1 1 1]; % b-spline curves

noPts      = 60;

hold on
nrbctrlplot(originalCurve);
axis equal


opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,'splinecurve.eps',opts)

% generating points on offset curve
thickness = 0.1;

xi=linspace(0,1,noPts);
offsetPts = zeros(noPts,2);

for i=1:length(xi)
    xii = xi(i);
    [N dNdxi] = NURBS1DBasisDers(xii,p,knots,weights);
    pts = controlPts(1:3,:);
    if xii > 0.5
        pts = controlPts(2:4,:);
    end
    x   = N    * pts;
    dx  = dNdxi* pts;
    dx  = dx/norm(dx);
    offsetPts(i,:) = x + [-dx(2) dx(1)]*thickness;
end

plot(offsetPts(:,1),offsetPts(:,2),'r*-','LineWidth',1.8);

% initial guess curve = a line

% cpoints(1,:) = offsetPts(1,:);
% cpoints(4,:) = offsetPts(end,:);
% 
% cpoints(2,:) = cpoints(1,:)+[1 0];
% cpoints(3,:) = cpoints(1,:)+[2 0];
% 
% cpoints0 = cpoints;


cpoints        = controlPts(:,1:2);

cpoints(1,  :) = offsetPts(1,:);
cpoints(end,:) = offsetPts(end,:);

vec = cpoints(end,:) - cpoints(1,:);

dlam = 1/(3);
for i=1:2   
    cpoints(1+i,:) = cpoints(1,:)+ dlam*i*vec;    
end

cpoints0 = cpoints;

samPts = computePoints(knots,weights,cpoints,xi);

plot(samPts(:,1),samPts(:,2),'bs-','LineWidth',1.8);

energy = getEnergy(offsetPts,samPts);

% disturbe the guess curve to compute numerically the gradient of the
% energy

h      = 1e-8;
eps    = 0.001;
eps1   = 0.01;
error  = 10;
error1 = 10;
dgamma = 0.0001;
gamma  = 0;
iiter  = 1;

while error > eps
    samPts   = computePoints(knots,weights,cpoints,xi);
    f        = getEnergy(offsetPts,samPts);
    
    cpoints(2,:) = cpoints(2,:) + [h 0];
    samPts   = computePoints(knots,weights,cpoints,xi);
    energy1x = getEnergy(offsetPts,samPts);
    cpoints(2,:) = cpoints(2,:) + [0 h] - [h 0];
    samPts   = computePoints(knots,weights,cpoints,xi);
    energy1y = getEnergy(offsetPts,samPts);
    
    dir = 1/h*[energy1x-energy energy1y-energy];
    %dir = dir/norm(dir);
    grad(1,:) = dir;
    
    cpoints(2,:) = cpoints0(2,:);
    cpoints(3,:) = cpoints(3,:) + [h 0];
    samPts   = computePoints(knots,weights,cpoints,xi);
    energy1x = getEnergy(offsetPts,samPts);
    cpoints(3,:) = cpoints(3,:) + [0 h] - [h 0];
    samPts   = computePoints(knots,weights,cpoints,xi);
    energy1y = getEnergy(offsetPts,samPts);
    
    dir = 1/h*[energy1x-energy energy1y-energy];
    %dir = dir;
    grad(2,:) = dir;
    
    s = 1;
    for k=1:10 
        cpnew = cpoints0;
        cpnew(2:3,:) = cpnew(2:3,:) - grad*s;        
        samPts   = computePoints(knots,weights,cpnew,xi);
        fnew     = getEnergy(offsetPts,samPts);
        if (fnew < f + s*alpha*(-grad)'*grad) 
            break;  
        else
            s = s*beta; 
        end
    end
    
    cpoints0(2:3,:) = cpoints0(2:3,:) - grad*s;
    cpoints = cpoints0;
    
    samPts   = computePoints(knots,weights,cpoints,xi);
    energy    = getEnergy(offsetPts,samPts);
    error = energy;
    
    disp (sprintf(' %s %i %s %5.4e ', 'Iter',iiter, ':', error) );
    iiter = iiter + 1;
end

figure
hold on
plot(offsetPts(:,1),offsetPts(:,2),'r--','LineWidth',1.3);

cp        = zeros(4,4);
cp(1:2,:) = cpoints0';
cp(4,:)   = 1;
offsetCurve = nrbmak(cp,knots);
nrbctrlplot(offsetCurve);

plot(cpoints0(:,1),cpoints0(:,2),'r-o',...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','g',...
    'MarkerSize',10);

