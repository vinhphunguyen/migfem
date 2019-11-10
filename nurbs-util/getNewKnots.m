function newKnots = getNewKnots(uKnot,p,weights,controlPts,dir,elConnU)

% Insert new knots into uKnot in such a way that the resulting mesh
% in the physical space is uniform.
% Vinh Phu Nguyen
% Delft University of Technology, March 2012

uniqueKnot = unique(uKnot);
noPts      = length(uniqueKnot);
node       = zeros(noPts,2);

eps = 1e-12;

%%
%% the following is problem-dependent!!!

if (dir==1)
    index = find(abs(controlPts(:,2)+6)<=eps);
else
    index = find(controlPts(:,1)==0);
end

points = controlPts(index,:);
%%

for i=1:noPts
    xi        = uniqueKnot(i);
    node(i,1) = NURBSinterpolation(xi,p, uKnot, points(:,1), weights);
    node(i,2) = NURBSinterpolation(xi,p, uKnot, points(:,2), weights);
end

newKnots = [];

for i=1:noPts-1
    if dir==1
        x1 = node(i,  1);
        x2 = node(i+1,1)
    else
        x1 = node(i,  2);
        x2 = node(i+1,2)
    end
    
    x    = 0.5*(x1 + x2);
    
    xi   = inverseMapping1DNURBS (x,uKnot,points,p,dir,elConnU,1);
    
    newKnots = [newKnots xi];
end
