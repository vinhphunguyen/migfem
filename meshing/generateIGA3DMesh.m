%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% MESH GENERATION: connectivity ...
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

uniqueUKnots   = unique(uKnot);
uniqueVKnots   = unique(vKnot);
uniqueWKnots   = unique(wKnot);

noElemsU       = length(uniqueUKnots)-1; % # of elements xi dir.
noElemsV       = length(uniqueVKnots)-1; % # of elements eta dir.
noElemsW       = length(uniqueWKnots)-1; % # of elements zeta dir.

%%%%%%%%%%%%%%%%%%%%%%%%%%
% chan =
%        1 2 3 4
%        5 6 7 8
%        9 10 11 12
%        13 14 15 16
% for a 4x2x2 control points


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

% combine info from two directions to build the elements
% element is numbered as follows
%  5 | 6 | 7 | 8
% ---------------
%  1 | 2 | 3 | 4
% for a 4x2 mesh

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






