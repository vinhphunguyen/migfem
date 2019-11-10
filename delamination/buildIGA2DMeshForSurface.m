function [mesh] = buildIGA2DMeshForSurface (surface)
%
% From the surface object created using nrbmak (...)
% convert to MIGFEM data structure

dim = surface.dim-1;


[p,q,uKnot,vKnot,noPtsX,noPtsY,weights,controlPts] = ...
               convertNURBSSurface (surface);
           
uniqueUKnots   = unique(uKnot);
uniqueVKnots   = unique(vKnot);
noElemsU       = length(uniqueUKnots)-1; % # of elements xi dir.
noElemsV       = length(uniqueVKnots)-1; % # of elements eta dir.
chan           = zeros(noPtsY,noPtsX);
count = 1;

for i=1:noPtsY
    for j=1:noPtsX
        chan(i,j) = count;
        count = count + 1;
    end
end

% determine our element ranges and the corresponding 
% knot indices along each direction

[elRangeU,elConnU] = buildConnectivity(p,uKnot,noElemsU);
[elRangeV,elConnV] = buildConnectivity(q,vKnot,noElemsV);
noElems = noElemsU * noElemsV;
element = zeros(noElems,(p+1)*(q+1));

e = 1;
for v=1:noElemsV
    vConn = elConnV(v,:);
    for u=1:noElemsU
        c = 1;
        uConn = elConnU(u,:);
        for i=1:length(vConn)
            for j=1:length(uConn)
              element(e,c) = chan(vConn(i),uConn(j));
              c = c + 1;
            end
        end
        e = e + 1;
    end  
end

index = zeros(noElems,2);
count = 1;

for j=1:size(elRangeV,1)
    for i=1:size(elRangeU,1)
        index(count,1) = i;
        index(count,2) = j;        
        count = count + 1;
    end
end



[C,Cxi,Cet]  = bezierExtraction2D(uKnot,vKnot,p,q);

mesh.p          = p;
mesh.q          = q;
mesh.uKnot      = uKnot;
mesh.vKnot      = vKnot;
mesh.noPtsX     = noPtsX;
mesh.noPtsY     = noPtsY;
mesh.noPts      = noPtsY*noPtsX;
mesh.weights    = weights;
mesh.controlPts = controlPts(:,1:dim);
mesh.C          = C;
mesh.locElems   = element;
mesh.globElems  = element;
mesh.noElemsU   = noElemsU;
mesh.noElemsV   = noElemsV;
mesh.noElems    = noElems;
mesh.index      = index;
mesh.elRangeU   = elRangeU;
mesh.elRangeV   = elRangeV;
mesh.elConnU    = elConnU;
mesh.elConnV    = elConnV;





