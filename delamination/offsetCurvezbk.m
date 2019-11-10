function [offsetCurve,offsetPts] = offsetCurve(originalCurve,thickness,...
                                               alpha,beta,eps,maxIter)

% computing the offset of a NURBS curve
% Inputs:
%  - originalCurve: the base curve
%  - thickness:     offset distance
%  - alpha,beta:    backtracking line search params
%  - eps:           tolerance
%  - maxIter:       after maxIter iterations, accept the solution anyway.
%
% Algorithm: gradient descent method with backtracking linesearch.
% Vinh Phu Nguyen, nvinhphu@gmail.com
% Cardiff University, Wales, UK, April 2013

%% get data from the structure originalCurve
knotVec    = originalCurve.knots;
p          = originalCurve.order-1;
controlPts = originalCurve.coefs(1:2,:)';
weights    = originalCurve.coefs(4,:);
noCtrPts   = size(controlPts,1);

%% computing offset points

% how to distribute these offset points?
% First option: uniform distribution
% noPts     = 200;
% xi        = linspace(0,1,noPts);
% Second option: divide knots into spans, for each span using an uniform
% distribution.

%xi=unique(xi); % remove duplicated values

offsetPts     = [];
xi            = [];
uKnotVec      = unique(knotVec);
noElems       = length(uKnotVec)-1;
multiplicity  = histc(knotVec,uKnotVec); % find duplicated knots
[ielements]   = buildGA1DMesh (knotVec,p);

figure
hold on

% loop over knot spans

noPts = [30 60 30];

specialPts = [];

for ik=1:noElems
    elemXi = linspace(uKnotVec(ik),uKnotVec(ik+1),noPts(ik))';
    %if ik ~= 1
    elemXi(end)=[];
    %end
    
    sctr   = ielements(ik,:);
    pnts   = [];
    for i=1:length(elemXi)
        xii = elemXi(i);
        
        [N dNdxi] = NURBS1DBasisDers(xii,p,knotVec,weights);
        s   = findspan(noCtrPts-1,p,xii,knotVec);
        pts = controlPts(s-1:s-1+p,:);
        %pts1 = controlPts(sctr,:); % not correct for evaluation
        %exactly at knots
        x   = N    * pts;
        dx  = dNdxi* pts;
        dx  = dx/norm(dx);
        pnt = x + [-dx(2) dx(1)]*thickness;
        pnts = [pnts;pnt];
    end
    
    id0 = [];
    
    % if xii is a knot in knotVec, then check its multiplicity
    xii = elemXi(1);
    ff = find(abs(uKnotVec-xii)<1e-10);
    
    if ~isempty(ff)
        mult = multiplicity(ff);
    end
    
    % knot with C^0 continuity
    if (mult == p)
        xii
        [N dNdxi] = NURBS1DBasisDers(xii-1e-4,p,knotVec,weights);
        pts  = controlPts(ielements(ik-1,:),:);
        x1   = N    * pts;
        dx1  = dNdxi* pts;
        dx1  = dx1/norm(dx1);
        
        n1   = [-dx1(2) dx1(1)];
        
        %            sctr   = ielements(ik+1,:);
        
        [N dNdxi] = NURBS1DBasisDers(xii,p,knotVec,weights);
        s   = findspan(noCtrPts-1,p,xii,knotVec);
        pts = controlPts(s-1:s-1+p,:);
        x2   = N    * pts;
        dx2  = dNdxi* pts;
        dx2  = dx2/norm(dx2);
        n2   = [-dx2(2) dx2(1)];
        averagedNormal = n1+n2;
        averagedNormal = averagedNormal/norm(averagedNormal);
        costheta = dot(n1,averagedNormal);
        pnt = x2 + averagedNormal*thickness/costheta;        
        specialPts = [specialPts; pnt];
        id  = find(pnts(:,1)<pnt(1));
        id0 = find(offsetPts(:,1)>pnt(1));
        pnts(id,:) = [];
        elemXi(id) = [];
        
    end
    
    if ~isempty(id0)
    offsetPts(id0,:) = [];
    xi(id0) = [];
    offsetPts(end,:) = pnt;
    xi(end)  = xii;
    end
    
    xi = [xi;elemXi];
    offsetPts = [offsetPts; pnts];
    
    
    mult = 100;
end

xi;

%xi        = linspace(0,1,length(xi));

%% gradient descent method
% initial guess curve = a line

cpoints        = controlPts(:,1:2);


% the first and last control points are the first and last offset points

cpoints(1,  :) = offsetPts(1,:);
cpoints(end,:) = offsetPts(end,:);

% if there are C^0 points, then control points there are known as well

cpoints(3,:)   = specialPts(1,:);
cpoints(5,:)   = specialPts(2,:);

fixedPntIds    = [1 3 5 7];
freePntIds     = setdiff(1:noCtrPts,fixedPntIds);

plot(specialPts(:,1),specialPts(:,2),'yd');

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

h      = 1e-10;
error  = 10;
iiter  = 1;

while error > eps
    samPts   = computePoints(knotVec,weights,cpoints,xi);
    f        = getEnergy(offsetPts,samPts);
    
    % determine the gradients numerically
    
    for ip = 1:length(freePntIds)
        i = freePntIds(ip);
        cpoints(i,:) = cpoints(i,:) + [h 0];
        samPts   = computePoints(knotVec,weights,cpoints,xi);
        energy1x = getEnergy(offsetPts,samPts);
        cpoints(i,:) = cpoints(i,:) - [h 0] + [0 h];
        samPts   = computePoints(knotVec,weights,cpoints,xi);
        energy1y = getEnergy(offsetPts,samPts);
        
        dir = 1/h*[energy1x-energy energy1y-energy];
        grad(ip,:) = dir;
        cpoints(i,:) = cpoints0(i,:);
    end
    
    s = 1;
    for k=1:10
        cpnew = cpoints0;
        cpnew(freePntIds,:) = cpnew(freePntIds,:) - grad*s;
        samPts   = computePoints(knotVec,weights,cpnew,xi);
        fnew     = getEnergy(offsetPts,samPts);
        if (fnew < f + s*alpha*(-grad)'*grad)
            break;
        else
            s = s*beta;
        end
    end
    
    cpoints0(freePntIds,:) = cpoints0(freePntIds,:) - grad*s;
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

% cpoints(4,2)=1.58;
% samPts   = computePoints(knotVec,weights,cpoints,xi);
% energy   = getEnergy(offsetPts,samPts)
    
%cpoints0(4,2)=1.58;
cp        = zeros(4,size(controlPts,1));
cp(1:2,:) = cpoints0';
cp(4,:)   = weights;
% cp(4,4) = 1.1;
% cp(1:2,4) = cp(1:2,4)*1.1;
offsetCurve = nrbmak(cp,knotVec);

%% debugging


plot(offsetPts(:,1),offsetPts(:,2),'r--','LineWidth',1.3);
plot(samPts(:,1),   samPts(:,2),'cy--','LineWidth',1.3);
size(samPts)
size(offsetPts)
nrbctrlplot(originalCurve);
axis equal
nrbctrlplot(offsetCurve);



