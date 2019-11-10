function ret = generateMesh(obj,node_pattern,node_patternLocal)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% MESH GENERATION: connectivity ...
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

uniqueUKnots   = unique(obj.uKnot);
uniqueVKnots   = unique(obj.vKnot);

noPtsX         = length(obj.uKnot)-obj.p-1;
noPtsY         = length(obj.vKnot)-obj.q-1;
noElemsU       = length(uniqueUKnots)-1; % # of elements xi dir.
noElemsV       = length(uniqueVKnots)-1; % # of elements eta dir.

% determine our element ranges and the corresponding 
% knot indices along each direction

[obj.elRangeU,obj.elConnU] = buildConnectivity(obj.p,obj.uKnot,noElemsU);
[obj.elRangeV,obj.elConnV] = buildConnectivity(obj.q,obj.vKnot,noElemsV);

% combine info from two directions to build the elements
% element is numbered as follows
%  5 | 6 | 7 | 8
% ---------------
%  1 | 2 | 3 | 4 
% for a 4x2 mesh

noElems = noElemsU * noElemsV;

element = zeros(noElems,(obj.p+1)*(obj.q+1));

e = 1;
for v=1:noElemsV
    vConn = obj.elConnV(v,:);
    for u=1:noElemsU
        c = 1;
        uConn = obj.elConnU(u,:);
        for i=1:length(vConn)
            for j=1:length(uConn)
              element(e,c) = node_pattern(vConn(i),uConn(j));
              c = c + 1;
            end
        end
        e = e + 1;
    end
end

elementLocal = zeros(noElems,(obj.p+1)*(obj.q+1));

e = 1;
for v=1:noElemsV
    vConn = obj.elConnV(v,:);
    for u=1:noElemsU
        c = 1;
        uConn = obj.elConnU(u,:);
        for i=1:length(vConn)
            for j=1:length(uConn)
              elementLocal(e,c) = node_patternLocal(vConn(i),uConn(j));
              c = c + 1;
            end
        end
        e = e + 1;
    end
end

index = zeros(noElems,2);
count = 1;

for j=1:size(obj.elRangeV,1)
    for i=1:size(obj.elRangeU,1)
        index(count,1) = i;
        index(count,2) = j;
        
        count = count + 1;
    end
end

obj.element      = element;
obj.elementLocal = elementLocal;
obj.index        = index;

end




