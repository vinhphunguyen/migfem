function xi = inverseMapping1DNURBS (x0,knots,controlPts,p,dir,elConnU,mode)
% inverse mapping of univariate NURBS interpolation.
% Input: 
%  x0        : point on the NURBS curve
%  knots     : knot vector
%  p         : basis order      
%  controlPts: control points define this curve
%  dir       : dir=1 for x and dir=2 for y
%  connU     : connectivity
% Method:
%  (1) find span in which x0 belongs to
%  (2) Newton-Raphson method to find xi


global weights

% knots      = [0 0 0 0.5 1 1 1];
% controlPts = [0 0; 1 2; 2 2; 3 0];
% p          = 2;
% 
% weights   = [1 1 1 1]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute knots images

uniqueKnot = unique(knots);
noPts      = length(uniqueKnot);
node       = zeros(noPts,2);

for i=1:noPts
  xi = uniqueKnot(i);  
  node(i,1) = NURBSinterpolation(xi,p, knots, controlPts(:,1), weights);
  node(i,2) = NURBSinterpolation(xi,p, knots, controlPts(:,2), weights);    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find in which element x0 belongs to
% x0;
% dir;
% node(:,dir);
node
tem    = node(:,dir) - x0;
if mode == 1
    elemId = length(find(tem<=0));
else
    elemId = length(find(tem>=0));
end

if elemId > length(elConnU)
    elemId = length(elConnU);
end

sctr   = elConnU(elemId,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Newton-Raphson iteration

noIterMax = 10;
epsilon   = 1e-4;
inc       = 1;
count     = 1;

% initial value of xi

xi        = uniqueKnot(elemId);

while (inc < noIterMax)
    [N dNdxi] = NURBS1DBasisDers(xi,p,knots,weights);        
    x         = N     * controlPts(sctr,dir);    
    df        = dNdxi * controlPts(sctr,dir);  
    f         = x - x0;
    xi        = xi - f/df;

    if (abs(f) < epsilon)
        inc  = noIterMax + 1;        
        %disp(['inverse mapping converged in ', num2str(count), ' iterations'])
        %disp(['with residual ', num2str(abs(f))])
    else
        inc = inc + 1;
    end    
    count = count + 1;
end
