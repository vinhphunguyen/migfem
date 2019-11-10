function [mesh] = buildIGA3DMesh(solid)

[p,q,r,uKnot,vKnot,wKnot,noPtsX,noPtsY,noPtsZ,weights,controlPts] = ...
               convert3DNURBS (solid);
           
uniqueUKnots   = unique(uKnot);
uniqueVKnots   = unique(vKnot);
uniqueWKnots   = unique(wKnot);

noElemsU       = length(uniqueUKnots)-1; % # of elements xi dir.
noElemsV       = length(uniqueVKnots)-1; % # of elements eta dir.
noElemsW       = length(uniqueWKnots)-1; % # of elements zeta dir.

chan  = zeros(noPtsZ,noPtsY,noPtsX);
count = 1;

for i=1:noPtsZ
    for j=1:noPtsY
        for k=1:noPtsX
            chan(i,j,k) = count;
            count       = count + 1;
        end
    end
end

% determine our element ranges and the corresponding 
% knot indices along each direction

[elRangeU,elConnU] = buildConnectivity(p,uKnot,noElemsU);
[elRangeV,elConnV] = buildConnectivity(q,vKnot,noElemsV);
[elRangeW,elConnW] = buildConnectivity(r,wKnot,noElemsW);

noElems = noElemsU * noElemsV * noElemsW;
element = zeros(noElems,(p+1)*(q+1)*(r+1));
e = 1;
for w=1:noElemsW
    wConn = elConnW(w,:);
    for v=1:noElemsV
        vConn = elConnV(v,:);
        for u=1:noElemsU
            c = 1;
            uConn = elConnU(u,:);
            
            for i=1:length(wConn)
                for j=1:length(vConn)
                    for k=1:length(uConn)
                        element(e,c) = chan(wConn(i),vConn(j),uConn(k));
                        c = c + 1;
                    end
                end
            end
            e = e + 1;
        end
    end
end

index = zeros(noElems,3);
count = 1;

for i=1:size(elRangeW,1)
    for j=1:size(elRangeV,1)
        for k=1:size(elRangeU,1)
            index(count,1) = k;
            index(count,2) = j;
            index(count,3) = i;
            count = count + 1;
        end
    end
end


[C,Cxi,Cet,Cze] = bezierExtraction3D(uKnot,vKnot,wKnot,p,q,r);
Cxe             = bezierExtraction2D(uKnot,vKnot,p,q);

mesh.p          = p;
mesh.q          = q;
mesh.r          = r;
mesh.uKnot      = uKnot;
mesh.vKnot      = vKnot;
mesh.wKnot      = wKnot;
mesh.noPts      = size(controlPts,1);
mesh.noPtsX     = noPtsX;
mesh.noPtsY     = noPtsY;
mesh.noPtsZ     = noPtsZ;
mesh.weights    = weights;
mesh.controlPts = controlPts(:,1:3);
mesh.C          = C;
mesh.Cxe        = Cxe;
mesh.locElems   = element;
mesh.globElems  = element;
mesh.rangeU     = elRangeU;
mesh.rangeV     = elRangeV;
mesh.rangeW     = elRangeW;
mesh.index      = index;
mesh.noElemsU   = noElemsU;
mesh.noElemsV   = noElemsV;
mesh.noElemsW   = noElemsW;
mesh.elConnU    = elConnU;
mesh.elConnV    = elConnV;
mesh.elConnW    = elConnW;
mesh.noElems    = noElems;


