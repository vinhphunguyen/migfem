function solid = doKRefinementSolid(solid,ods,refCount)

pnew = ods(1);
qnew = ods(2);
knew = ods(3);

refCountX = refCount(1);
refCountY = refCount(2);
refCountZ = refCount(3);

orders = solid.order - 1;
p      = orders(1);
q      = orders(2);
k      = orders(3);

solid  = nrbdegelev(solid,[pnew-p qnew-q knew-k]); 

uKnot  = cell2mat(solid.knots(1));
vKnot  = cell2mat(solid.knots(2));
wKnot  = cell2mat(solid.knots(3));

for i=1:refCountX
    uKnotVectorU = unique(uKnot);          
    newKnotsX = uKnotVectorU(1:end-1) + 0.5*diff(uKnotVectorU);   
    newKnots  = {newKnotsX [] []};    
    % h-refinement    
    solid     = nrbkntins(solid,newKnots);    
    uKnot     = cell2mat(solid.knots(1));
end

for i=1:refCountY
    vKnotVectorU = unique(vKnot);         
    newKnotsY = vKnotVectorU(1:end-1) + 0.5*diff(vKnotVectorU);   
    newKnots  = {[] newKnotsY []};    
    % h-refinement    
    solid     = nrbkntins(solid,newKnots);    
    vKnot     = cell2mat(solid.knots(2));
end

for i=1:refCountZ
    wKnotVectorU = unique(wKnot);          
    newKnotsZ = wKnotVectorU(1:end-1) + 0.5*diff(wKnotVectorU);   
    newKnots  = {[] [] newKnotsZ};    
    % h-refinement    
    solid     = nrbkntins(solid,newKnots);    
    wKnot     = cell2mat(solid.knots(3));
end
