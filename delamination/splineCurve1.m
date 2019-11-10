addpath ../C_files/
addpath ../fem_util/
addpath ../examples/

clear all

knots      = [0 0 0 1 1 1];

uniqueKnot = unique(knots);

controlPts  = [0 1; 1 1; 1 0];
weights     = [1 1/sqrt(2) 1]; % b-spline curves

cp = zeros(4,3);
cp(1:2,:) = controlPts';
cp(4,:)   = 1;
cp(4,2)   = 1/sqrt(2);
cp(1:2,2) = cp(1:2,2) * 1/sqrt(2);

originalCurve = nrbmak(cp,knots);

p          = 2;
noPts      = 80;

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
    x   = N    * pts;
    dx  = dNdxi* pts;
    dx  = dx/norm(dx);
    offsetPts(i,:) = x + [-dx(2) dx(1)]*thickness;
end

plot(offsetPts(:,1),offsetPts(:,2),'r*-','LineWidth',1.8);

% initial guess curve = a line

cpoints(1,:) = offsetPts(1,:);
cpoints(3,:) = offsetPts(end,:);

cpoints(2,:) = 0.5*(cpoints(1,:) + cpoints(3,:)); %controlPts(2,:);


cpoints0 = cpoints;

samPts = computePoints(knots,weights,cpoints,xi);

plot(samPts(:,1),samPts(:,2),'bs-','LineWidth',1.8);

energy = getEnergy(offsetPts,samPts);

% disturbe the guess curve to compute numerically the gradient of the
% energy

h      = 1e-8;
eps    = 0.00001;
eps1   = 0.00001;
error  = 10;
error1 = 10;
dgamma = 0.0005;
gamma  = 0;
iiter  = 1;

while error > eps
    cpoints(2,:) = cpoints(2,:) + [h 0];
    samPts   = computePoints(knots,weights,cpoints,xi);
    energy1x = getEnergy(offsetPts,samPts);
    cpoints(2,:) = cpoints(2,:) + [0 h] - [h 0];
    samPts   = computePoints(knots,weights,cpoints,xi);
    energy1y = getEnergy(offsetPts,samPts);
    
    dir = 1/h*[energy1x-energy energy1y-energy];
    dir = dir/norm(dir);
    grad(1,1:2) = dir;
    
    while error > eps1
        gamma = gamma + dgamma;        
        cpoints0(2,:) = cpoints0(2,:) - grad*dgamma;
        samPts   = computePoints(knots,weights,cpoints0,xi);
        error    = getEnergy(offsetPts,samPts);
    end
    
    cpoints = cpoints0;
    
    disp (sprintf(' %s %i %s %5.4e ', 'Iter',iiter, ':', error) );
    iiter = iiter + 1;
    
    cp        = zeros(4,3);
    cp(1:2,:) = cpoints0';
    cp(4,:)   = 1;
    originalCurve = nrbmak(cp,knots);
    nrbctrlplot(originalCurve);
end

figure
hold on
plot(offsetPts(:,1),offsetPts(:,2),'r-','LineWidth',1.8);

cp        = zeros(4,3);
cp(1:2,:) = cpoints';
cp(4,:)   = weights;
cp(1,:) = cp(1,:) .* weights;
cp(2,:) = cp(2,:) .* weights;
originalCurve = nrbmak(cp,knots);
nrbctrlplot(originalCurve);

% plot(cpoints0(:,1),cpoints0(:,2),'r-o',...
%     'MarkerEdgeColor','k',...
%     'MarkerFaceColor','g',...
%     'MarkerSize',10);

