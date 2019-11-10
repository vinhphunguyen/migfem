% h-refinement using knot insertion
% knots are inserted in such a way that a uniform mesh
% in physical space is obtained. Positions of new knots
% are determined by calling function "getNewKnots".
% Vinh Phu Nguyen, March 2012
% Delft University of Technology
% Adapted from ISOGAT, A.V. Vuong.

for c=1:refineCount
    uKnotVectorU = unique(uKnot);
    uKnotVectorV = unique(vKnot);
    
    noElemsU       = length(uKnotVectorU)-1; % # of elements xi dir.
    noElemsV       = length(uKnotVectorV)-1; % # of elements eta dir.
    
    [elRangeU,elConnU] = buildConnectivity(p,uKnot,noElemsU);
    [elRangeV,elConnV] = buildConnectivity(q,vKnot,noElemsV);

    % new knots along two directions
    
    newKnotsX = getNewKnots(uKnot,p,weights,controlPts,1,elConnU);
    newKnotsY = getNewKnots(vKnot,q,weights,controlPts,2,elConnV);
    
    %% h-refinement (NURBS) in x-direction
    dim = size(controlPts,2);
    
    nonewkX      = size(newKnotsX,2);
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
    
    %% h-refinement (NURBS) in y-direction)
    nonewkY       = size(newKnotsY,2);
    newprojcoord  = zeros(noPtsX*noPtsY+nonewkY*noPtsY,dim+1);
       
    for i=1:noPtsX
        %create index for reading controlPoints
        rcpindex     = i:noPtsX:noPtsX*noPtsY;
        locCP        = controlPts(rcpindex,:);
        locweights   = weights(rcpindex);
        locprojcoord = nurb2proj(noPtsY, locCP, locweights);
        
        % refinement of y
        [tempknotVectorY,tempcontrolPoints] = ...
            RefineKnotVectCurve(noPtsY-1,q,vKnot,locprojcoord,newKnotsY,nonewkY-1);
        
        wcpindex = i:noPtsX:noPtsX*(noPtsY+nonewkY);
        newprojcoord(wcpindex,:) = tempcontrolPoints;
    end
    
    vKnot                 = tempknotVectorY;
    [controlPts, weights] = proj2nurbs(newprojcoord);
    noPtsY                = noPtsY + nonewkY;
end








