%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% MESH GENERATION: connectivity ...
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Vinh Phu Nguyen, March 2012
% nvinhphu@gmail.com

uniqueUKnots   = unique(uKnot);
uniqueVKnots   = unique(vKnot);

noElemsU       = length(uniqueUKnots)-1; % # of elements xi dir.
noElemsV       = length(uniqueVKnots)-1; % # of elements eta dir.

%%%%%%%%%%%%%%%%%%%%%%%%%%
% chan = 
%        1 2 3 4
%        5 6 7 8
%        9 10 11 12
% for a 4x3 control points


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

% combine info from two directions to build the elements
% element is numbered as follows
%  5 | 6 | 7 | 8
% ---------------
%  1 | 2 | 3 | 4 
% for a 4x2 mesh

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

% ptsIndex = zeros(noElems,(p+1)*(q+1));
% 
% for e=1:noElems
%    idu    = index(e,1);
%    idv    = index(e,2);
%    connU  = elConnU(idu,:);
%    connV  = elConnV(idv,:);
%    
%    noFnsU = length(connU);
%    noFnsV = length(connV);
%    
%    c = 1;
%    for i=1:noFnsU
%        for j=1:noFnsV
%            ptsIndex(e,c) = (connV(j)-1)*noPtsX + connU(i);
%            c = c + 1;
%        end
%    end
% end




