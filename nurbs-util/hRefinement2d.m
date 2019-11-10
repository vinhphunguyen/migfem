% clear all
% clc
% 
% % Data
% 
% %plateHoleData  % plate with a hole
% %plate1ElemData % rectangular plate in tension
% %plate2ElemData
% annularData
% %LShapedData
% %LShapedC1Data
% 
% 
% 
% refineCount = 2;

% Vinh Phu Nguyen
% Delft University of Technology
% Adapted from ISOGAT, A.V. Vuong.

for c=1:refineCount
    uKnotVectorU = unique(uKnot);
    uKnotVectorV = unique(vKnot);
    
%     if ( c == 1)
%         noPtsX         = size(controlPts,1);
%         noPtsY         = size(controlPts,2);
%         controlPts     = reshape(controlPts,noPtsX*noPtsY,2);
%         weights        = ones(1,noPtsX*noPtsY)';
%     end
    
    
    % new knots along two directions
    
    newKnotsX = uKnotVectorU(1:end-1) + 0.5*diff(uKnotVectorU);
    newKnotsY = uKnotVectorV(1:end-1) + 0.5*diff(uKnotVectorV);
    
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








