% original data

p = 2;
q = 1;


uKnot = [0 0 0 1 1 1];
vKnot = [0 0 1 1];

controlPts=[-1 0;-1 1;0 1;
    -2 0;-2 2; 0 2];

noPtsX = 3;
noPtsY = 2;

weights = ones(1,noPtsX*noPtsY)';

weights(2)=1/sqrt(2);
weights(5)=1/sqrt(2);

% add new knots

newKnotsX = [0.5 0.5];


%% h-refinement (NURBS) in x-direction
dim = size(controlPts,2);

nonewkX      = size(newKnotsX,2)
newprojcoord = zeros(noPtsX*noPtsY+nonewkX*noPtsY,dim+1);

rstart = 1;
wstart = 1;

for j=1:noPtsY
    rstop = rstart + noPtsX-1;
    wstop = wstart + noPtsX-1 + nonewkX;
    
    locCP        = controlPts(rstart:rstop,:);
    locweights   = weights   (rstart:rstop);
    locprojcoord = nurb2proj(noPtsX, locCP, locweights);
    
    % refinement of x
    [tempknotVectorX,tempControlPts] = ...
        RefineKnotVectCurve(noPtsX-1,p,uKnot,locprojcoord,newKnotsX,nonewkX-1);
    
    newprojcoord(wstart:wstop,:)=tempControlPts;
    wstart = wstop+1;
    rstart = rstop+1;
end

uKnot                 = tempknotVectorX;
[controlPts, weights] = proj2nurbs(newprojcoord);
noPtsX                = noPtsX+nonewkX;




