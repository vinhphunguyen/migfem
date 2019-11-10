function [offsetCurve,offsetPts] = offsetCurveNormal(originalCurve,thickness,alpha,beta,eps,maxIter)

% computing the offset of a NURBS curve
% Inputs:
%  - originalCurve: the base curve 
%  - thickness:     offset distance
%  - alpha,beta:    backtracking line search params
%  - eps:           tolerance
%
% Algorithm: gradient descent method
% Vinh Phu Nguyen, nvinhphu@gmail.com
% Cardiff University, Wales, UK, April 2013

knotVec    = originalCurve.knots;
p          = originalCurve.order-1; 
controlPts = originalCurve.coefs(1:2,:)';
weights    = originalCurve.coefs(4,:);
noCtrPts   = size(controlPts,1);

% our controlPts only stores (x,y,z) not (w*x,w*y,w*z)

controlPts(:,1) = controlPts(:,1)./weights';
controlPts(:,2) = controlPts(:,2)./weights';

% computing offset points

% how to distribute these offset points?
% First option: uniform distribution 
% noPts     = 200;
% xi        = linspace(0,1,noPts);
% Second option: divide knots into spans, for each span using an uniform
% distribution.

xi = [];
uKnotVec = unique(knotVec);

for ik=1:length(uKnotVec)-1
    xi = [xi; linspace(uKnotVec(ik),uKnotVec(ik+1),40)'];
end

xi=unique(xi);

offsetPts = zeros(length(xi),2);

for i=1:length(xi)
    xii = xi(i);
    [N dNdxi] = NURBS1DBasisDers(xii,p,knotVec,weights);    
    s   = findspan(noCtrPts-1,p,xii,knotVec);
    pts = controlPts(s-p+1:s+1,:);
    x   = N    * pts;
    dx  = dNdxi* pts;
    dx  = dx/norm(dx);
    offsetPts(i,:) = x + [-dx(2) dx(1)]*thickness;
end

%% gradient descent method 
% initial guess curve = a line

cpoints        = controlPts(:,1:2);

cpoints(1,  :) = offsetPts(1,:);
cpoints(end,:) = offsetPts(end,:);

vec = cpoints(end,:) - cpoints(1,:);

% dlam = 1/(noCtrPts-1);
% for i=1:noCtrPts-2    
%     cpoints(1+i,:) = cpoints(1,:)+ dlam*i*vec;    
% end

cpoints0 = cpoints;
samPts = computePoints(knotVec,weights,cpoints,xi);
energy = getEnergy(offsetPts,samPts);

% disturbe the guess curve to compute numerically the gradient of the
% energy

h      = 1e-8;
error  = 10;
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
    
    if iiter > maxIter 
        break 
    end
end

cp        = zeros(4,size(controlPts,1));
cp(1:2,:) = cpoints0';
cp(4,:)   = weights
cp(1,:) = cp(1,:) .* weights;
cp(2,:) = cp(2,:) .* weights;
offsetCurve = nrbmak(cp,knotVec);

%% debugging 


figure
hold on
plot(offsetPts(:,1),offsetPts(:,2),'r--','LineWidth',1.3);
%plot(samPts(:,1),   samPts(:,2),'cy--','LineWidth',1.3);

nrbctrlplot(originalCurve);
axis equal
nrbctrlplot(offsetCurve);



