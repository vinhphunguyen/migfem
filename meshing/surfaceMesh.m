function [element,index] = surfaceMesh (uKnot,vKnot,surfaceNodes,p,q,...
                           noPtsX,noPtsY,elRangeU,elRangeV,elConnU,elConnV)

% Build a NURBS surface mesh from the control points defining
% this surface. This is useful for imposing boundary conditions
% on a surface of a NURBS solid.
% It is based on the file generate2DMesh.m 
% Vinh Phu Nguyen
% Delft University of Technology
                       
uniqueUKnots   = unique(uKnot);
uniqueVKnots   = unique(vKnot);

noElemsU       = length(uniqueUKnots)-1; % # of elements xi dir.
noElemsV       = length(uniqueVKnots)-1; % # of elements eta dir.

%%%%%%%%%%%%%%%%%%%%%%%%%%
% chan = 

chan = reshape(surfaceNodes,noPtsX,noPtsY)';

%chan = surfaceNodes;

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






